package com.compomics.util.experiment.io.identification.idfilereaders;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.aminoacids.AminoAcid;
import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidSequence;
import com.compomics.util.experiment.biology.atoms.Atom;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.mass_spectrometry.SpectrumProvider;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.io.IoUtil;
import com.compomics.util.io.flat.SimpleFileReader;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.waiting.WaitingHandler;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.stream.Collectors;
import javax.xml.bind.JAXBException;
import org.xmlpull.v1.XmlPullParser;
import org.xmlpull.v1.XmlPullParserException;
import org.xmlpull.v1.XmlPullParserFactory;

/**
 * Simple IdfileReader for Pepxml files.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class PepxmlIdfileReader implements IdfileReader {

    /**
     * List of the spectrum matches in the file.
     */
    private ArrayList<SpectrumMatch> spectrumMatches = null;
    /**
     * The name of the search engine which was used to create the file.
     */
    private String searchEngine = null;
    /**
     * The version of the search engine which was used to create the file.
     */
    private String searchEngineVersion = null;
    /**
     * The file to parse.
     */
    private File idFile;
    /**
     * The name of the spectrum file.
     */
    private String spectrumFileName;
    /**
     * Stores the mass differences of the fixed modifications. The key is the
     * amino acid residue as a single upper case character and the element is
     * the list of the mass differences of the masses targeting that residue.
     */
    private HashMap<Character, ArrayList<Double>> fixedModificationsMassDiffs;
    /**
     * Stores the masses of the fixed modifications.
     */
    private ArrayList<Double> fixedModificationMasses;
    /**
     * Stores the masses of the fixed n-terminal modifications.
     */
    private ArrayList<Double> fixedNTerminalModifications = new ArrayList<>();
    /**
     * Stores the masses of the fixed c-terminal modifications.
     */
    private ArrayList<Double> fixedCTerminalModifications = new ArrayList<>();

    /**
     * Blank constructor for instantiation purposes.
     */
    public PepxmlIdfileReader() {

    }

    /**
     * Constructor.
     *
     * @param idFile the file to parse
     */
    public PepxmlIdfileReader(
            File idFile
    ) {
        this.idFile = idFile;
    }

    /**
     * Parses the identification file.
     *
     * @param waitingHandler A waiting handler returning information about the
     * progress and allowing canceling the parsing.
     * @param expandAaCombinations If true the combinations of amino acids will
     * be expanded.
     * @param overwriteExtension If true, the extension of the input file will
     * be overwritten to mgf.
     * @param spectrumProvider A spectrum provider with the spectra of the file loaded.
     *
     * @throws XmlPullParserException Exception thrown if an error occurred while parsing the xml file.
     * @throws IOException Exception thrown if an error occurred while reading the file.
     */
    private void parseFile(
            WaitingHandler waitingHandler, 
            boolean expandAaCombinations, 
            boolean overwriteExtension,
            SpectrumProvider spectrumProvider
    ) throws XmlPullParserException, IOException {

        // Create the pull parser.
        XmlPullParserFactory factory = XmlPullParserFactory.newInstance(System.getProperty(XmlPullParserFactory.PROPERTY_NAME), null);
        factory.setNamespaceAware(true);
        XmlPullParser parser = factory.newPullParser();

        // Create a reader for the input file.
        try ( SimpleFileReader reader = SimpleFileReader.getFileReader(idFile)) {
            
            // Set the XML Pull Parser to read from this reader.
            parser.setInput(reader.getReader());
            
            // Start the parsing.
            int type;
            boolean hasMatch = false;

            HashMap<Long, SpectrumMatch> spectrumMatchesMap = new HashMap<>();
            spectrumMatches = new ArrayList<>();
            SpectrumMatch currentMatch = null;
            Integer currentCharge = null;

            // Go through the whole document.
            while ((type = parser.next()) != XmlPullParser.END_DOCUMENT) {
                String tagName = parser.getName();
                if (type == XmlPullParser.START_TAG && tagName.equals("msms_run_summary")) {
                    
                    parseRunSummary(parser, overwriteExtension);
                    
                }
                if (type == XmlPullParser.START_TAG && tagName.equals("search_summary")) {
                    parseSearchSummary(parser);
                }
                if (type == XmlPullParser.START_TAG && tagName.equals("spectrum_query")) {
                    
                    currentMatch = parseSpectrumQuery(
                            parser, 
                            spectrumProvider
                    );
                    SpectrumMatch previousMatch = spectrumMatchesMap.get(currentMatch.getKey());
                    if (previousMatch != null) {
                        currentMatch = previousMatch;
                    }
                    for (int i = 0; i < parser.getAttributeCount(); i++) {
                        String attributeName = parser.getAttributeName(i);
                        if (attributeName.equals("assumed_charge")) {
                            String value = parser.getAttributeValue(i);
                            try {
                                currentCharge = new Integer(value.trim());
                            } catch (Exception e) {
                                throw new IllegalArgumentException("Charge " + value + " could not be parsed. Integer expected.");
                            }
                        }
                    }
                }
                if (type == XmlPullParser.START_TAG && tagName.equals("search_hit")) {
                    if (currentMatch == null) {
                        throw new IllegalArgumentException("No spectrum match when parsing search hit.");
                    }
                    if (currentCharge == null) {
                        throw new IllegalArgumentException("No charge found when parsing search hit of spectrum " + currentMatch.getKey() + ".");
                    }
                    PeptideAssumption peptideAssumption = parseSearchHit(parser, currentCharge);
                    Peptide peptide = peptideAssumption.getPeptide();
                    String peptideSequence = peptide.getSequence();
                    hasMatch = true;
                    boolean found = false;
                    if (currentMatch.getAllPeptideAssumptions().count() > 0) {
                        for (PeptideAssumption tempPeptideAssumption : currentMatch.getAllPeptideAssumptions().collect(Collectors.toList())) {
                            Peptide tempPeptide = tempPeptideAssumption.getPeptide();
                            if (peptide.getSequence().equals(tempPeptide.getSequence())) {
                                boolean sameModifications = peptide.getNVariableModifications() == tempPeptide.getNVariableModifications();
                                if (sameModifications) {
                                    for (ModificationMatch originalMatch : peptide.getVariableModifications()) {
                                        boolean modFound = false;
                                        for (ModificationMatch otherMatch : tempPeptide.getVariableModifications()) {
                                            if (originalMatch.getModification().equals(otherMatch.getModification()) && originalMatch.getSite() == otherMatch.getSite()) {
                                                modFound = true;
                                                break;
                                            }
                                        }
                                        if (!modFound) {
                                            sameModifications = false;
                                            break;
                                        }
                                    }
                                }
                                if (sameModifications) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (!found) {

                        Advocate advocate = Advocate.getAdvocate(searchEngine);
                        if (expandAaCombinations && AminoAcidSequence.hasCombination(peptideSequence)) {

                            ModificationMatch[] previousModificationMatches = peptide.getVariableModifications();

                            for (StringBuilder expandedSequence : AminoAcidSequence.getCombinations(peptide.getSequence())) {

                                ModificationMatch[] newModificationMatches = Arrays.stream(previousModificationMatches)
                                        .map(modificationMatch -> modificationMatch.clone())
                                        .toArray(ModificationMatch[]::new);

                                Peptide newPeptide = new Peptide(expandedSequence.toString(), newModificationMatches, true);
                                PeptideAssumption newAssumption = new PeptideAssumption(newPeptide, peptideAssumption.getRank(),
                                        peptideAssumption.getAdvocate(), peptideAssumption.getIdentificationCharge(),
                                        peptideAssumption.getScore(), peptideAssumption.getIdentificationFile());
                                currentMatch.addPeptideAssumption(advocate.getIndex(), newAssumption);
                            }
                        } else {
                            currentMatch.addPeptideAssumption(advocate.getIndex(), peptideAssumption);
                        }
                    }
                }

                if (type == XmlPullParser.END_TAG && tagName.equals("spectrum_query")) {

                    if (hasMatch) {

                        long key = currentMatch.getKey();

                        if (!spectrumMatchesMap.containsKey(key)) {

                            spectrumMatchesMap.put(key, currentMatch);
                            spectrumMatches.add(currentMatch);

                        }

                        hasMatch = false;
                        currentMatch = null;
                        currentCharge = null;

                    }
                }
                
                if (waitingHandler != null) {
                    waitingHandler.isRunCanceled();
                    
                }
            }

            spectrumMatchesMap.clear();

        }
    }

    /**
     * Parses a search hit.
     *
     * @param parser the XML parser
     * @param charge the charge of the hit
     *
     * @return the peptide assumption in the search hit
     *
     * @throws XmlPullParserException Exception thrown if an error occurred while parsing the xml file.
     * @throws IOException Exception thrown if an error occurred while reading the file.
     */
    private PeptideAssumption parseSearchHit(
            XmlPullParser parser, 
            int charge
    ) throws XmlPullParserException, IOException {

        Integer rank = null;
        String sequence = null;
        ArrayList<ModificationMatch> modificationMatches = new ArrayList<>();
        Double score = null;

        for (int i = 0; i < parser.getAttributeCount(); i++) {
            String name = parser.getAttributeName(i);
            if (name.equals("hit_rank")) {
                String value = parser.getAttributeValue(i);
                try {
                    rank = new Integer(value.trim());
                } catch (Exception e) {
                    throw new IllegalArgumentException("An error occurred while parsing rank " + value + ". Integer expected.");
                }
            } else if (name.equals("peptide")) {
                sequence = parser.getAttributeValue(i).trim();
            }
        }

        int type;
        while ((type = parser.next()) != XmlPullParser.START_TAG) {
        }

        String tagName = parser.getName();
        if (tagName.equals("modification_info")) {

            for (int i = 0; i < parser.getAttributeCount(); i++) {
                String attributeName = parser.getAttributeName(i);
                if (attributeName.equals("mod_nterm_mass") || attributeName.equals("mod_cterm_mass")) {

                    String value = parser.getAttributeValue(i).trim();
                    Double terminalMass = null;
                    try {
                        terminalMass = new Double(value);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("An error occurred while parsing modification terminal mass " + value + ". Number expected.");
                    }

                    // check if the terminal modification is fixed or variable
                    boolean variableModification;
                    if (attributeName.equals("mod_nterm_mass")) {
                        variableModification = !fixedNTerminalModifications.contains(terminalMass);
                    } else {
                        variableModification = !fixedCTerminalModifications.contains(terminalMass);
                    }

                    if (variableModification) {

                        int site;
                        if (attributeName.equals("mod_nterm_mass")) {
                            site = 1;
                            terminalMass -= Atom.H.getMonoisotopicMass();
                        } else { // c-term
                            site = sequence.length();
                            terminalMass -= (Atom.O.getMonoisotopicMass() + Atom.H.getMonoisotopicMass());

                            // fix for older comet pepxml files
                            if (searchEngine != null && searchEngine.equalsIgnoreCase("Comet")
                                    && searchEngineVersion != null
                                    && !searchEngineVersion.equalsIgnoreCase("2015.02 rev. 4")
                                    && !searchEngineVersion.equalsIgnoreCase("2015.02 rev. 5")) { // @TODO: make more generic...
                                terminalMass -= Atom.H.getMonoisotopicMass();
                            }
                        }

                        char aa = sequence.charAt(site - 1);
                        terminalMass = Util.roundDouble(terminalMass, 2);
                        String tempModificationName = terminalMass + "@" + aa;
                        ModificationMatch modificationMatch = new ModificationMatch(tempModificationName, site);
                        modificationMatches.add(modificationMatch);

                    }
                }
            }

            while ((type = parser.next()) != XmlPullParser.END_DOCUMENT) {
                tagName = parser.getName();
                if (tagName != null) {
                    if (tagName.equals("mod_aminoacid_mass")) {
                        Integer site = null;
                        for (int i = 0; i < parser.getAttributeCount(); i++) {
                            String attributeName = parser.getAttributeName(i);
                            if (attributeName.equals("position")) {
                                String value = parser.getAttributeValue(i);
                                try {
                                    site = new Integer(value);
                                } catch (Exception e) {
                                    throw new IllegalArgumentException("An error occurred while parsing modification position " + value + ". Integer expected.");
                                }
                            }
                        }
                        if (site != null) {
                            Double modifiedAaMass = null;
                            for (int i = 0; i < parser.getAttributeCount(); i++) {
                                String attributeName = parser.getAttributeName(i);
                                if (attributeName.equals("mass")) {
                                    String value = parser.getAttributeValue(i);
                                    try {
                                        modifiedAaMass = new Double(value);
                                    } catch (Exception e) {
                                        throw new IllegalArgumentException("An error occurred while parsing modification mass " + value + ". Number expected.");
                                    }
                                }
                            }

                            if (modifiedAaMass != null) {
                                char aa = sequence.charAt(site - 1);
                                AminoAcid aminoAcid = AminoAcid.getAminoAcid(aa);

                                // see if the amino acid also has a fixed modification
                                //
                                // example:
                                //  carbamidomethyl _and_ pyrolidone from carbamidomethylated c:
                                //
                                //  <modification_info modified_peptide="C[143]EQALLQVAK">
                                //      <mod_aminoacid_mass position="1" mass="143.004100"/>
                                //  </modification_info>
                                //
                                double fixedModificationMass = 0;
                                boolean variableModification;

                                if (fixedModificationMasses.contains(modifiedAaMass)) {
                                    variableModification = false;
                                } else {
                                    if (fixedModificationsMassDiffs.get(aa) != null) {
                                        for (Double tempMassDiff : fixedModificationsMassDiffs.get(aa)) {
                                            fixedModificationMass += tempMassDiff;
                                        }
                                    }
                                    variableModification = true;
                                }

                                if (variableModification) {
                                    double modificationMass = modifiedAaMass - fixedModificationMass - aminoAcid.getMonoisotopicMass();
                                    modificationMass = Util.roundDouble(modificationMass, 2);
                                    String tempModificationName = modificationMass + "@" + aa;
                                    ModificationMatch modificationMatch = new ModificationMatch(tempModificationName, site);
                                    modificationMatches.add(modificationMatch);
                                }
                            }
                        }
                    } else if (type == XmlPullParser.END_TAG && parser.getName().equals("modification_info")) {
                        while ((type = parser.next()) != XmlPullParser.START_TAG) {
                        }
                        break;
                    }
                }
            }
        }

        while (type != XmlPullParser.END_DOCUMENT) {
            tagName = parser.getName();
            if (tagName != null) {
                if (type == XmlPullParser.START_TAG && parser.getName().equals("search_score")) {
                    String name = null;
                    String value = null;
                    for (int i = 0; i < parser.getAttributeCount(); i++) {
                        String attributeName = parser.getAttributeName(i);
                        if (attributeName.equals("name")) {
                            name = parser.getAttributeValue(i);
                        } else if (attributeName.equals("value")) {
                            value = parser.getAttributeValue(i);
                        }
                    }

                    if (name != null && value != null) {
                        if (name.equals("expect") || name.equals("Morpheus Score")) { // @TODO: what about the IdentiPy hyperscore (IdentiPy:hyperscore (MS:1002989)) and/or IdentiPy:RHNS (MS:1002988)?
                            try {
                                score = new Double(value);
                            } catch (Exception e) {
                                throw new IllegalArgumentException("Impossible to parse expectation value " + value + ". Number expected.");
                            }
                        }
                    }
                } else if (type == XmlPullParser.END_TAG && tagName.equals("search_hit")) {
                    break;
                }
            }
            type = parser.next();
        }

        Peptide peptide = new Peptide(
                sequence, 
                modificationMatches.toArray(new ModificationMatch[modificationMatches.size()]), 
                true
        );
        Advocate advocate = Advocate.getAdvocate(searchEngine);
        
        return new PeptideAssumption(
                peptide, 
                rank, 
                advocate.getIndex(), 
                charge, 
                score, 
                idFile.getName()
        );
    }

    /**
     * Parses a spectrum query.
     *
     * @param parser The XML parser.
     * @param spectrumProvider A spectrum provider with the spectra of the file loaded.
     *
     * @return the spectrum match in this spectrum query
     *
     * @throws XmlPullParserException Exception thrown if an error occurred while parsing the xml file.
     * @throws IOException Exception thrown if an error occurred while reading the file.
     */
    private SpectrumMatch parseSpectrumQuery(
            XmlPullParser parser,
            SpectrumProvider spectrumProvider
    ) throws XmlPullParserException, IOException {

        Integer index = null;
        String spectrumId = null;
        String spectrumNativeID = null;

        for (int i = 0; i < parser.getAttributeCount(); i++) {
            String name = parser.getAttributeName(i);
            if (name.equals("spectrum")) {
                spectrumId = parser.getAttributeValue(i);
            } else if (name.equals("index")) {
                String value = parser.getAttributeValue(i);
                try {
                    index = new Integer(value.trim());
                } catch (Exception e) {
                    throw new IllegalArgumentException("An error occurred while parsing index " + value + ". Integer expected.");
                }
            } else if (name.equals("spectrumNativeID")) {
                spectrumNativeID = parser.getAttributeValue(i);
            }
        }

        if (index == null) {
            throw new IllegalArgumentException("No index found for spectrum " + spectrumId + ".");
        }

        String spectrumTitle;

        if (spectrumNativeID != null) {
            
            spectrumTitle = spectrumNativeID;
        
        } else {
            
            spectrumTitle = spectrumProvider.getSpectrumTitles(spectrumFileName)[index];
            
        }
        
        SpectrumMatch spectrumMatch = new SpectrumMatch(spectrumFileName, spectrumTitle);

        return spectrumMatch;
    }

    /**
     * Parses the run summary.
     *
     * @param parser the XML parser
     * @param overwriteExtension if true, the extension of the input file will
     * be overwritten to mgf
     *
     * @throws XmlPullParserException Exception thrown if an error occurred while parsing the xml file.
     * @throws IOException Exception thrown if an error occurred while reading the file.
     */
    private void parseRunSummary(
            XmlPullParser parser, 
            boolean overwriteExtension
    ) throws XmlPullParserException, IOException {

        // Something like  <msms_run_summary base_name="D:\path\filename" raw_data="extention"> is expected 
        String path = "";
        for (int i = 0; i < parser.getAttributeCount(); i++) {
            String name = parser.getAttributeName(i);
            if (name.equals("base_name")) {
                path += parser.getAttributeValue(i);
            } else if (!overwriteExtension && name.equals("raw_data")) {
                path += parser.getAttributeValue(i);
            }
        }

        if (overwriteExtension) {
            path += ".mgf";
        }

        File spectrumFile = new File(path);
        spectrumFileName = IoUtil.getFileName(spectrumFile);
        
    }

    /**
     * Parses the search summary.
     *
     * @param parser the XML parser
     *
     * @throws XmlPullParserException Exception thrown if an error occurred while parsing the xml file.
     * @throws IOException Exception thrown if an error occurred while reading the file.
     */
    private void parseSearchSummary(
            XmlPullParser parser
    ) throws XmlPullParserException, IOException {
        
        for (int i = 0; i < parser.getAttributeCount(); i++) {
            String name = parser.getAttributeName(i);
            if (name.equals("search_engine")) {
                searchEngine = parser.getAttributeValue(i);
            } else if (name.equals("search_engine_version")) {
                searchEngineVersion = parser.getAttributeValue(i);
            }
        }

        // extract the required information about the modifications
        fixedModificationsMassDiffs = new HashMap<>();
        fixedModificationMasses = new ArrayList<>();
        fixedNTerminalModifications = new ArrayList<>();
        fixedCTerminalModifications = new ArrayList<>();

        int type;

        while ((type = parser.next()) != XmlPullParser.END_DOCUMENT) {

            if (type == XmlPullParser.END_TAG && parser.getName() != null) {
                if (parser.getName().equals("search_summary")) {
                    break;
                }
            }

            if (type == XmlPullParser.START_TAG) {

                String tagName = parser.getName();

                if (type == XmlPullParser.START_TAG && tagName.equals("aminoacid_modification")) {

                    Character aminoacid = null;
                    Boolean variable = null;
                    Double massDiff = null;
                    Double mass = null;

                    for (int i = 0; i < parser.getAttributeCount(); i++) {
                        String name = parser.getAttributeName(i);
                        if (name.equals("aminoacid")) {
                            aminoacid = parser.getAttributeValue(i).charAt(0);
                        } else if (name.equals("massdiff")) {
                            massDiff = new Double(parser.getAttributeValue(i));
                        } else if (name.equals("mass")) {
                            mass = new Double(parser.getAttributeValue(i));
                        } else if (name.equals("variable")) {
                            String variableAsString = parser.getAttributeValue(i);
                            if (variableAsString.equalsIgnoreCase("Y")) {
                                variable = true;
                            } else if (variableAsString.equalsIgnoreCase("N")) {
                                variable = false;
                            }
                        }
                    }

                    if (variable != null && massDiff != null && mass != null && aminoacid != null) {
                        if (!variable) {
                            ArrayList<Double> massDiffs = fixedModificationsMassDiffs.get(aminoacid);
                            if (massDiffs == null) {
                                massDiffs = new ArrayList<>();
                            }
                            massDiffs.add(massDiff);
                            fixedModificationsMassDiffs.put(aminoacid, massDiffs);
                            fixedModificationMasses.add(mass);
                        }
                    } else {
                        throw new IllegalArgumentException("An error occurred while parsing aminoacid_modification element. Missing values.");
                    }

                } else if (type == XmlPullParser.START_TAG && tagName.equals("terminal_modification")) {

                    Boolean variable = null;
                    Double mass = null;
                    String terminus = null;

                    for (int i = 0; i < parser.getAttributeCount(); i++) {
                        String name = parser.getAttributeName(i);
                        if (name.equals("terminus")) {
                            String terminusAsString = parser.getAttributeValue(i);
                            if (terminusAsString.equalsIgnoreCase("N") || terminusAsString.equalsIgnoreCase("C")) {
                                terminus = terminusAsString;
                            }
                        } else if (name.equals("mass")) {
                            mass = new Double(parser.getAttributeValue(i));
                        } else if (name.equals("variable")) {
                            String variableAsString = parser.getAttributeValue(i);
                            if (variableAsString.equalsIgnoreCase("Y")) {
                                variable = true;
                            } else if (variableAsString.equalsIgnoreCase("N")) {
                                variable = false;
                            }
                        }
                    }

                    if (variable != null && mass != null && terminus != null) {
                        if (!variable) {
                            if (terminus.equalsIgnoreCase("N")) {
                                fixedNTerminalModifications.add(mass);
                            } else {
                                fixedCTerminalModifications.add(mass);
                            }
                        }
                    } else {
                        throw new IllegalArgumentException("An error occurred while parsing terminal_modification element. Missing values.");
                    }
                }
            }
        }
    }

    @Override
    public HashMap<String, ArrayList<String>> getSoftwareVersions() {
        ArrayList<String> versions = new ArrayList<>(1);
        versions.add(searchEngineVersion);
        HashMap<String, ArrayList<String>> result = new HashMap<>(1);
        result.put(searchEngine, versions);
        return result;
    }

    @Override
    public String getExtension() {
        return ".pep.xml";
    }

    @Override
    public void close() throws IOException {
        // nothing to do here
    }

    @Override
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler, 
            SearchParameters searchParameters
    ) 
            throws IOException, SQLException, ClassNotFoundException, InterruptedException, JAXBException, XmlPullParserException {
        
        return getAllSpectrumMatches(
                spectrumProvider, 
                waitingHandler, 
                searchParameters, 
                null, 
                true
        );
    }

    @Override
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler, 
            SearchParameters searchParameters,
            SequenceMatchingParameters sequenceMatchingPreferences, 
            boolean expandAaCombinations
    ) 
            throws IOException, SQLException, ClassNotFoundException, InterruptedException, JAXBException, XmlPullParserException {
        
        if (spectrumMatches == null) {
        
            parseFile(
                    waitingHandler, 
                    expandAaCombinations, 
                    true, 
                    spectrumProvider
            );
        }
        
        return spectrumMatches;
    }

    @Override
    public boolean hasDeNovoTags() {
        return false;
    }
}
