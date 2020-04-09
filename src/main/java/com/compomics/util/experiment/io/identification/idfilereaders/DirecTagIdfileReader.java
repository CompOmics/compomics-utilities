package com.compomics.util.experiment.io.identification.idfilereaders;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.aminoacids.AminoAcid;
import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidSequence;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.experiment.identification.spectrum_assumptions.TagAssumption;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.amino_acid_tags.Tag;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
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
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Set;
import javax.xml.bind.JAXBException;

/**
 * An identification file reader for Direct tag results.
 *
 * @author Marc Vaudel
 */
public class DirecTagIdfileReader implements IdfileReader {

    /**
     * The name of the tags generator used to create the file.
     */
    private String tagsGenerator;
    /**
     * The version of the tags generator.
     */
    private String tagsGeneratorVersion;
    /**
     * The copyright.
     */
    private String copyRight;
    /**
     * The license.
     */
    private String license;
    /**
     * The time of sequencing start.
     */
    private String timeStart;
    /**
     * The time of sequencing end.
     */
    private String timeEnd;
    /**
     * The tagging time.
     */
    private Double taggingTimeSeconds;
    /**
     * The number of processing nodes.
     */
    private Integer nProcessingNode;
    /**
     * The file used as input.
     */
    private String inputFile;
    /**
     * The tags parameters in a map.
     */
    private final HashMap<String, String> tagsParameters = new HashMap<>();
    /**
     * Returns the content of the columns for a spectrum line. Name &gt; index
     * in the column.
     */
    private final HashMap<String, Integer> spectrumLineContent = new HashMap<>();
    /**
     * Returns the content of the columns for a tag line. Name &gt; index in the
     * column.
     */
    private final HashMap<String, Integer> tagLineContent = new HashMap<>();
    /**
     * The file inspected.
     */
    private File tagFile;
    /**
     * The mass to add to the C-terminal gap so that is corresponds to a peptide
     * fragment.
     */
    public final double cTermCorrection = 0;
    /**
     * The mass to add to the N-terminal gap so that is corresponds to a peptide
     * fragment.
     */
    public final double nTermCorrection = 0;
    /**
     * The DirecTag parameters.
     */
    private DirecTagParameters direcTagParameters;
    /**
     * The residues modified by the different PTMs in the DynamicMods tag.
     * Indexed on the symbol used to represent the PTM.
     */
    private HashMap<Character, Character> dynamicModsResidues;

    /**
     * Default constructor for the purpose of instantiation.
     */
    public DirecTagIdfileReader() {

    }

    /**
     * Constructors, parses a file.
     *
     * @param tagFile the file to parse
     */
    public DirecTagIdfileReader(
            File tagFile
    ) {

        this.tagFile = tagFile;

    }

    /**
     * Returns the name of the different parameters names found.
     *
     * @return the name of the different parameters names found
     */
    public Set<String> getTagsParametersNames() {
        return tagsParameters.keySet();
    }

    /**
     * Returns the tagging parameter corresponding to a given parameter name.
     *
     * @param tagParameterName the name of the parameter of interest
     *
     * @return the parameter of interest
     */
    public String getTagParameter(
            String tagParameterName
    ) {

        return tagsParameters.get(tagParameterName);

    }

    /**
     * Parses the parameters section.
     *
     * @param reader The file reader.
     *
     * @return true if the end of the file was reached
     *
     * @throws IOException if an IOException occurs
     */
    private boolean parseParameters(
            SimpleFileReader reader
    ) {

        String line;
        while ((line = reader.readLine()) != null) {

            if (line == null || line.startsWith("H	TagsParameters")) {

                break;

            } else if (line == null) {

                throw new IllegalArgumentException("Unexpected end of file while parsing the parameters.");

            } else if (line.startsWith("H(S)") || line.startsWith("H(T)") || line.startsWith("S") || line.startsWith("T")) {

                throw new IllegalArgumentException("Unexpected end of parameters section.");

            } else {

                line = line.substring(1).trim();

                if (line.startsWith("TagsGeneratorVersion")) {

                    tagsGeneratorVersion = line.substring(line.indexOf("\t")).trim();

                } else if (line.startsWith("TagsGenerator")) {

                    tagsGenerator = line.substring(line.indexOf("\t")).trim();

                } else if (line.contains("(c)")) {

                    copyRight = line;

                } else if (line.contains("License")) {

                    license = line;

                } else if (line.startsWith("Tagging started at")) {

                    tagsGeneratorVersion = line.substring(line.indexOf("Tagging started at")).trim();

                } else if (line.startsWith("Tagging started at")) {

                    timeStart = line.substring(line.indexOf("Tagging started at")).trim();

                } else if (line.startsWith("Tagging finished at")) {

                    timeEnd = line.substring(line.indexOf("Tagging finished at")).trim();

                } else if (line.startsWith("Total tagging time:")) {

                    line = line.substring(line.indexOf(":") + 1).trim();
                    line = line.substring(0, line.indexOf(" ")).trim();

                    try {

                        taggingTimeSeconds = new Double(line);

                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else if (line.contains("node")) {

                    line = line.substring(line.indexOf(" ")).trim();
                    line = line.substring(0, line.indexOf(" ")).trim();

                    try {

                        nProcessingNode = new Integer(line);

                    } catch (Exception e) {
                        // ignore
                    }

                } else if (line.startsWith("InputFile")) {

                    inputFile = line.substring(line.indexOf("\t")).trim();

                }
            }
        }

        return line == null;

    }

    /**
     * Parses the tag parameters.
     *
     * @return true if the end of the file was reached
     */
    private boolean parseTagParameters(
            SimpleFileReader reader
    ) {

        String line;
        while ((line = reader.readLine()) != null) {

            if (line.trim().isEmpty()) {

                break;

            } else if (line == null) {

                throw new IllegalArgumentException("Unexpected end of file while parsing the tag parameters.");

            } else if (line.startsWith("H(S)") || line.startsWith("H(T)") || line.startsWith("S") || line.startsWith("T")) {

                throw new IllegalArgumentException("Unexpected end of tag parameters section.");

            } else {

                line = line.substring(1).trim();
                String[] components = line.split(", ");

                for (String component : components) {

                    int index = component.indexOf(": ");

                    if (index != -1) {

                        String key = component.substring(0, index).trim();
                        String value = component.substring(index + 1).trim();
                        tagsParameters.put(key, value);

                    }
                }
            }
        }

        return line == null;

    }

    /**
     * Parses the tables headers.
     *
     * @return true if the end of the file was reached
     */
    private boolean parseHeaders(
            SimpleFileReader reader
    ) {

        String line = reader.readLine();

        if (line != null) {

            parseHeaderLine(line);

        }

        line = reader.readLine();

        if (line != null) {

            parseHeaderLine(line);

        }

        return line == null;

    }

    /**
     * Parses a line corresponding to a header.
     *
     * @param linea line corresponding to a header
     */
    private void parseHeaderLine(
            String line
    ) {

        if (line.startsWith("S") || line.startsWith("T")) {

            throw new IllegalArgumentException("No Header found.");

        }

        if (line.startsWith("H(S)")) {

            line = line.substring(4).trim();
            String[] components = line.split("\t");

            for (int i = 0; i < components.length; i++) {

                spectrumLineContent.put(components[i], i);

            }

        } else if (line.startsWith("H(T)")) {

            line = line.substring(4).trim();
            String[] components = line.split("\t");

            for (int i = 0; i < components.length; i++) {

                tagLineContent.put(components[i], i);

            }
        }
    }

    /**
     * Sets the dynamic modifications from the tags parameters.
     */
    private void setDynamicMods() {

        // get the ptm residues from the DynamicMods field
        dynamicModsResidues = new HashMap<>();
        String dynamicMods = tagsParameters.get("DynamicMods"); // assume something like: "M 0 15.994915 N 1 0.984016 Q 2 0.984016"
        dynamicMods = dynamicMods.trim();
        if (!dynamicMods.isEmpty()) {
            String[] modElements = dynamicMods.split(" ");
            int index = 0;
            while (index + 2 < modElements.length) {
                dynamicModsResidues.put(modElements[index + 1].charAt(0), modElements[index].charAt(0));
                index += 3;
            }
        }
    }

    /**
     * Parses the results section.
     */
    private ArrayList<SpectrumMatch> parseResults(
            SimpleFileReader reader,
            String[] spectrumTitles
    ) {

        ArrayList<SpectrumMatch> result = new ArrayList<>();

        String spectrumFileName = IoUtil.getFileName(getInputFile());

        int spectrumCount = 0;
        Integer sIdColumnIndex = spectrumLineContent.get("ID");
        Integer chargeColumnIndex = spectrumLineContent.get("Charge");

        int lastId = -1, lastCharge = -1;
        int rank = 0;
        SpectrumMatch currentMatch = null;
        String line;

        while ((line = reader.readLine()) != null) {

            if (line.startsWith("S")) {

                int sId = ++spectrumCount;
                rank = 0;

                if (sIdColumnIndex != null) {

                    line = line.substring(1).trim();
                    String[] components = line.split("\t");
                    String id = components[sIdColumnIndex];
                    sId = Integer.parseInt(id.substring(id.indexOf("=") + 1));
                    String chargeString = components[chargeColumnIndex];
                    lastCharge = Integer.parseInt(chargeString);

                }
                if (sId != lastId) {

                    if (currentMatch != null && currentMatch.getAllTagAssumptions().count() > 0) {

                        result.add(currentMatch);
                    }

                    String spectrumTitle = spectrumTitles[sId];
                    currentMatch = new SpectrumMatch(spectrumFileName, spectrumTitle);
                    lastId = sId;

                }

            } else if (line.startsWith("T")) {

                ++rank;
                TagAssumption tagAssumption = getAssumptionFromLine(line, rank);
                //@TODO: check with the developers if this is correct
                tagAssumption.setIdentificationCharge(lastCharge);
                currentMatch.addTagAssumption(Advocate.direcTag.getIndex(), tagAssumption);

            }
        }

        if (currentMatch != null && currentMatch.getAllTagAssumptions().count() > 0) {

            result.add(currentMatch);
        }

        return result;

    }

    @Override
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler,
            SearchParameters searchParameters
    ) throws IOException, SQLException, ClassNotFoundException, InterruptedException, JAXBException {

        return getAllSpectrumMatches(
                spectrumProvider,
                waitingHandler,
                searchParameters,
                null,
                false
        );
    }

    @Override
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler,
            SearchParameters searchParameters,
            SequenceMatchingParameters sequenceMatchingPreferences,
            boolean expandAaCombinations
    ) throws IOException, IllegalArgumentException, SQLException, ClassNotFoundException, InterruptedException, JAXBException {

        direcTagParameters = (DirecTagParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.direcTag.getIndex());

        ArrayList<SpectrumMatch> result = new ArrayList<>(0);

        try ( SimpleFileReader reader = SimpleFileReader.getFileReader(tagFile)) {

            boolean endOfFile = parseParameters(reader);

            if (!endOfFile) {

                endOfFile = parseTagParameters(reader);

            }
            if (!endOfFile) {

                endOfFile = parseHeaders(reader);

            }
            if (!endOfFile) {

                setDynamicMods();
                result = parseResults(
                        reader,
                        spectrumProvider.getSpectrumTitles(
                                getInputFile().getName()
                        )
                );
            }
        }
        return result;
    }

    /**
     * Returns the assumption associated to a tag line. If a modification index
     * is found, an "X" is put in the tag sequence and a modification match
     * named after the given index is added.
     *
     * @param line the line
     * @param rank the rank of the assumption
     *
     * @return the assumption associated to a tag line
     */
    private TagAssumption getAssumptionFromLine(
            String line,
            int rank
    ) {

        line = line.substring(1).trim();
        String[] components = line.split("\t");
        Integer cGapIndex = tagLineContent.get("cTerminusMass");
        if (cGapIndex == null) {
            throw new IllegalArgumentException("Column cTerminusMass not found.");
        }
        Double cGap = new Double(components[cGapIndex]);
        if (cGap > 0 && cGap < cTermCorrection) {
            throw new IllegalArgumentException("Incompatible c-term gap " + cGap);
        } else if (cGap > 0) {
            cGap += cTermCorrection;
        }
        Integer nGapIndex = tagLineContent.get("nTerminusMass");
        if (nGapIndex == null) {
            throw new IllegalArgumentException("Column nTerminusMass not found.");
        }
        Double nGap = new Double(components[nGapIndex]);
        Integer tagIndex = tagLineContent.get("Tag");
        if (tagIndex == null) {
            throw new IllegalArgumentException("Column Tag not found.");
        }
        String tagSequence = components[tagIndex];
        StringBuilder residues = new StringBuilder(tagSequence.length());
        HashMap<Integer, ModificationMatch> modificationMatches = new HashMap<>();
        for (int i = 0; i < tagSequence.length(); i++) {
            char charAtI = tagSequence.charAt(i);
            try {
                AminoAcid aa = AminoAcid.getAminoAcid(charAtI);
                residues.append(aa.singleLetterCode);
            } catch (IllegalArgumentException e) {
                try {
                    // modified residue
                    String modIndexString = charAtI + "";
                    int modIndex = new Integer(modIndexString);
                    String utilitiesPtm = direcTagParameters.getUtilitiesModificationName(modIndex);
                    modificationMatches.put(i + 1, new ModificationMatch(utilitiesPtm, i + 1));
                    residues.append(dynamicModsResidues.get(modIndexString.charAt(0)));
                } catch (Exception e1) {
                    throw new IllegalArgumentException("No amino acid or modification could be mapped to tag component \"" + charAtI + "\" in tag \"" + tagSequence + "\".");
                }
            }
        }

        AminoAcidSequence tagAaSequence = new AminoAcidSequence(residues.toString());
        for (int i : modificationMatches.keySet()) {
            tagAaSequence.addVariableModification(modificationMatches.get(i));
        }
        Tag tag = new Tag(nGap, tagAaSequence, cGap);

        Integer chargeIndex = tagLineContent.get("TagChargeState");
        if (chargeIndex == null) {
            throw new IllegalArgumentException("Column TagChargeState not found.");
        }
        int charge = new Integer(components[chargeIndex]);

        Integer eValueIndex = tagLineContent.get("Total");
        if (eValueIndex == null) {
            throw new IllegalArgumentException("Column Total not found.");
        }
        double eValue = new Double(components[eValueIndex]);

        return new TagAssumption(Advocate.direcTag.getIndex(), rank, tag, charge, eValue);
    }

    /**
     * Returns the tags generator used to create the file.
     *
     * @return the tags generator used to create the file
     */
    public String getTagsGenerator() {
        return tagsGenerator;
    }

    /**
     * Returns the version of the tags generator used to create the file.
     *
     * @return the version of the tags generator used to create the file
     */
    public String getTagsGeneratorVersion() {
        return tagsGeneratorVersion;
    }

    /**
     * Returns the copyright.
     *
     * @return the copyright
     */
    public String getCopyRight() {
        return copyRight;
    }

    /**
     * Returns the license information of this file.
     *
     * @return the license information of this file
     */
    public String getLicense() {
        return license;
    }

    /**
     * Returns the starting time of the tagging as given in the file.
     *
     * @return the starting time of the tagging
     */
    public String getTimeStart() {
        return timeStart;
    }

    /**
     * Returns the ending time of the tagging as given in the file.
     *
     * @return the ending time of the tagging
     */
    public String getTimeEnd() {
        return timeEnd;
    }

    /**
     * Returns the tagging time in seconds as listed in the file.
     *
     * @return the tagging time in seconds as listed in the file
     */
    public Double getTaggingTimeSeconds() {
        return taggingTimeSeconds;
    }

    /**
     * Returns the number of processing nodes used.
     *
     * @return the number of processing nodes used
     */
    public Integer getnProcessingNode() {
        return nProcessingNode;
    }

    /**
     * Returns the spectrum file name as found in the parameters section.
     *
     * @return the spectrum file name
     */
    public File getInputFile() {
        return new File(inputFile);
    }

    @Override
    public String getExtension() {
        return ".tags";
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public HashMap<String, ArrayList<String>> getSoftwareVersions() {

        HashMap<String, ArrayList<String>> result = new HashMap<>();
        ArrayList<String> versions = new ArrayList<>();
        versions.add(tagsGeneratorVersion);
        result.put(tagsGenerator, versions);
        return result;
    }

    @Override
    public boolean hasDeNovoTags() {
        return true;
    }
}
