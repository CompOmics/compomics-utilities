package com.compomics.util.experiment.io.identifications.idfilereaders;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.AminoAcidPattern;
import com.compomics.util.experiment.biology.Atom;
import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.TagAssumption;
import com.compomics.util.experiment.identification.identification_parameters.PepnovoParameters;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.tags.Tag;
import com.compomics.util.experiment.io.identifications.IdfileReader;
import com.compomics.util.experiment.massspectrometry.Charge;
import com.compomics.util.experiment.massspectrometry.Spectrum;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.experiment.refinementparameters.PepnovoAssumptionDetails;
import com.compomics.util.waiting.WaitingHandler;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import uk.ac.ebi.pride.tools.braf.BufferedRandomAccessFile;

/**
 * This class can be used to parse PepNovo identification files.
 *
 * @author Marc Vaudel
 */
public class PepNovoIdfileReader extends ExperimentObject implements IdfileReader {

    /**
     * A map of all spectrum titles and the associated index in the random
     * access file.
     */
    private HashMap<String, Long> index;
    /**
     * The result file in random access.
     */
    private BufferedRandomAccessFile bufferedRandomAccessFile = null;
    /**
     * The name of the result file.
     */
    private String fileName;
    /**
     * The standard format.
     */
    public static final String tableHeader = "#Index	RnkScr	PnvScr	N-Gap	C-Gap	[M+H]	Charge	Sequence";

    /**
     * Default constructor for the purpose of instantiation.
     */
    public PepNovoIdfileReader() {

    }

    /**
     * Constructor, initiate the parser. Displays the progress using the waiting
     * handler. The close() method shall be used when the file reader is no
     * longer used.
     *
     * @param identificationFile the identification file to parse
     *
     * @throws FileNotFoundException exception thrown whenever the provided file
     * was not found
     * @throws IOException exception thrown whenever an error occurred while
     * reading the file
     */
    public PepNovoIdfileReader(File identificationFile) throws FileNotFoundException, IOException {
        this(identificationFile, null);
    }

    /**
     * Constructor, initiate the parser. Displays the progress using the waiting
     * handler. The close() method shall be used when the file reader is no
     * longer used.
     *
     * @param identificationFile the identification file to parse
     * @param waitingHandler a waiting handler providing progress feedback to
     * the user
     *
     * @throws FileNotFoundException exception thrown whenever the provided file
     * was not found
     * @throws IOException exception thrown whenever an error occurred while
     * reading the file
     */
    public PepNovoIdfileReader(File identificationFile, WaitingHandler waitingHandler) throws FileNotFoundException, IOException {

        bufferedRandomAccessFile = new BufferedRandomAccessFile(identificationFile, "r", 1024 * 100);
        fileName = Util.getFileName(identificationFile);

        if (waitingHandler != null) {
            waitingHandler.resetSecondaryProgressCounter();
            waitingHandler.setMaxSecondaryProgressCounter(100);
        }

        long progressUnit = bufferedRandomAccessFile.length() / 100;

        if (progressUnit == 0) {
            progressUnit = 1;
        }

        index = new HashMap<String, Long>();

        String line;
        while ((line = bufferedRandomAccessFile.readLine()) != null) {
            if (line.startsWith(">>")) {
                long currentIndex = bufferedRandomAccessFile.getFilePointer();

                String[] temp = line.split("\\s+");
                String formatted = "";
                for (int i = 3; i < temp.length; i++) {
                    formatted += (temp[i] + " ");
                }
                int endIndex = formatted.lastIndexOf("#Problem");
                if (endIndex == -1) {
                    endIndex = formatted.lastIndexOf("(SQS");
                }

                // Condition: Skip problematic spectra not containing (SQS) at the end of the line.
                if (endIndex > -1) {
                    String spectrumTitle = formatted.substring(0, endIndex).trim();
                    index.put(spectrumTitle, currentIndex);
                }
                if (waitingHandler != null) {
                    if (waitingHandler.isRunCanceled()) {
                        break;
                    }
                    waitingHandler.setSecondaryProgressCounter((int) (currentIndex / progressUnit));
                }

            }
        }
    }

    @Override
    public HashSet<SpectrumMatch> getAllSpectrumMatches(WaitingHandler waitingHandler) throws IOException, IllegalArgumentException, Exception {

        if (bufferedRandomAccessFile == null) {
            throw new IllegalStateException("The identification file was not set. Please use the appropriate constructor.");
        }

        HashSet<SpectrumMatch> spectrumMatches = new HashSet<SpectrumMatch>();

        if (waitingHandler != null) {
            waitingHandler.setSecondaryProgressCounterIndeterminate(false);
            waitingHandler.resetSecondaryProgressCounter();
            waitingHandler.setMaxSecondaryProgressCounter(index.size());
        }

        for (String title : index.keySet()) {

            String decodedTitle = URLDecoder.decode(title, "utf-8");
            SpectrumMatch currentMatch = new SpectrumMatch(Spectrum.getSpectrumKey(getMgfFileName(), decodedTitle));

            int cpt = 1;
            bufferedRandomAccessFile.seek(index.get(title));
            String line = bufferedRandomAccessFile.getNextLine().trim();
            boolean solutionsFound = true;
            if (line.startsWith("# No") || line.startsWith("# Charge") || line.startsWith("#Problem") || line.startsWith("# too")) {
                solutionsFound = false;
            } else if (!line.equals(tableHeader)) {
                throw new IllegalArgumentException("Unrecognized table format. Expected: \"" + tableHeader + "\", found:\"" + line + "\".");
            }

            while ((line = bufferedRandomAccessFile.getNextLine()) != null
                    && !line.equals("") && !line.startsWith(">>")) {
                currentMatch.addHit(Advocate.pepnovo.getIndex(), getAssumptionFromLine(line, cpt), true);
                cpt++;
            }
            if (solutionsFound) {
                spectrumMatches.add(currentMatch);
            }

            if (waitingHandler != null) {
                if (waitingHandler.isRunCanceled()) {
                    break;
                }
                waitingHandler.increaseSecondaryProgressCounter();
            }
        }

        return spectrumMatches;
    }

    /**
     * Returns the spectrum file name. This method assumes that the PepNovo
     * output file is the mgf file name + ".out"
     *
     * @return the spectrum file name
     */
    public String getMgfFileName() {
        return fileName.substring(0, fileName.length() - 4);
    }

    @Override
    public String getExtension() {
        return ".out";
    }

    @Override
    public void close() throws IOException {
        bufferedRandomAccessFile.close();
    }

    /**
     * Returns a Peptide Assumption from a PepNovo result line. The rank score
     * is taken as reference score. All additional parameters are attached as
     * PeptideAssumptionDetails.
     *
     * Note: Fixed PTMs are not annotated, variable PTMs are marked with the
     * PepNovo PTM tag (see PepnovoParameters to retrieve utilities names).
     *
     * @param line the line to parse
     * @param rank the rank of the assumption
     * @return the corresponding assumption
     */
    private TagAssumption getAssumptionFromLine(String line, int rank) {

        String[] lineComponents = line.trim().split("\t");

        Double rankScore = new Double(lineComponents[1]);
        Double pepNovoScore = new Double(lineComponents[2]);
        Double nGap = new Double(lineComponents[3]);
        Double cGap = new Double(lineComponents[4]);
        double cTermCorrection = Atom.O.getMonoisotopicMass() + Atom.H.getMonoisotopicMass() + 2 * ElementaryIon.proton.getTheoreticMass();
        if (cGap > 0 && cGap < cTermCorrection) {
            throw new IllegalArgumentException("Incompatible c-term gap " + cGap);
        } else if (cGap > 0) {
            cGap -= cTermCorrection;
        }
        Double mH = new Double(lineComponents[5]);
        Integer charge = new Integer(lineComponents[6]);
        String pepNovoSequence = lineComponents[7];
        String sequence = "";
        ArrayList<ModificationMatch> modificationMatches = new ArrayList<ModificationMatch>();
        String modificationMass = "", currentAA = "";
        int currentPtmLocation = 0;

        boolean nTermPtm = false;
        boolean cTermPtm = false;

        String ptmTag = "";

        // find and add the variable ptms
        for (int i = 0; i < pepNovoSequence.length(); i++) {
            String aa = pepNovoSequence.charAt(i) + "";

            if (aa.equals("^") || aa.equals("$")) {

                ptmTag = aa;

                if (aa.equals("^")) {
                    nTermPtm = true;
                } else {
                    cTermPtm = true;
                }

            } else {

                if (aa.equals("+") || aa.equals("-")) {
                    modificationMass += aa;
                } else {
                    try {
                        new Integer(aa);
                        modificationMass += aa;
                    } catch (Exception e) {
                        if (!modificationMass.equals("")) {

                            String pepNovoPtmTag = "";

                            if (nTermPtm || cTermPtm) {
                                pepNovoPtmTag += ptmTag;
                            } else {
                                pepNovoPtmTag += currentAA;
                            }

                            pepNovoPtmTag += modificationMass;
                            ModificationMatch modMatch = new ModificationMatch(pepNovoPtmTag, true, currentPtmLocation);
                            modMatch.setConfident(true);
                            modificationMatches.add(modMatch);
                            modificationMass = "";
                            nTermPtm = false;
                        }
                        AminoAcid aminoAcid = AminoAcid.getAminoAcid(aa);
                        if (aminoAcid == null) {
                            throw new IllegalArgumentException("Attempting to parse " + aa + " as amino acid in " + pepNovoSequence + ".");
                        }
                        sequence += aa;
                        currentAA = aa;
                        currentPtmLocation++;
                    }
                }
            }
        }

        if (!modificationMass.equals("")) {

            String pepNovoPtmTag = "";

            if (nTermPtm || cTermPtm) {
                pepNovoPtmTag += ptmTag;
            } else {
                pepNovoPtmTag += currentAA;
            }

            pepNovoPtmTag += modificationMass;

            ModificationMatch modMatch = new ModificationMatch(pepNovoPtmTag, true, currentPtmLocation);
            modificationMatches.add(modMatch);
        }

        AminoAcidPattern aminoAcidPattern = new AminoAcidPattern(sequence);
        for (ModificationMatch modificationMatch : modificationMatches) {
            aminoAcidPattern.addModificationMatch(modificationMatch.getModificationSite(), modificationMatch);
        }
        Tag tag = new Tag(nGap, aminoAcidPattern, cGap);
        TagAssumption tagAssumption = new TagAssumption(Advocate.pepnovo.getIndex(), rank, tag, new Charge(Charge.PLUS, charge), pepNovoScore);
        PepnovoAssumptionDetails pepnovoAssumptionDetails = new PepnovoAssumptionDetails();
        pepnovoAssumptionDetails.setRankScore(rankScore);
        pepnovoAssumptionDetails.setMH(mH);
        tagAssumption.addUrParam(pepnovoAssumptionDetails);

        return tagAssumption;
    }

    /**
     * Get a PTM.
     *
     * @param pepnovoParameters the pep novo parameters
     * @param pepNovoModification the PepNovo modification
     *
     * @return the PTM as a string
     */
    public static String getPTM(PepnovoParameters pepnovoParameters, String pepNovoModification) {

        Map<String, String> invertedPtmMap = pepnovoParameters.getPepNovoPtmMap();

        if (invertedPtmMap == null) {
            // @TODO: possible to rescue these?
            throw new IllegalArgumentException("Unsupported de novo search result. Please reprocess the data.");
        }

        String utilitesPtmName = invertedPtmMap.get(pepNovoModification);

        if (utilitesPtmName != null) {
            return utilitesPtmName;
        } else {
            throw new IllegalArgumentException("An error occurred while parsing the modification " + pepNovoModification + ".");
        }
    }

    @Override
    public String getSoftwareVersion() {
        throw new UnsupportedOperationException("Not supported yet."); //@TODO
    }
}
