package com.compomics.util.parameters.identification.search;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.tool_specific.AndromedaParameters;
import com.compomics.util.parameters.identification.tool_specific.CometParameters;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
import com.compomics.util.parameters.identification.tool_specific.MsAmandaParameters;
import com.compomics.util.parameters.identification.tool_specific.MsgfParameters;
import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.parameters.identification.tool_specific.NovorParameters;
import com.compomics.util.parameters.identification.tool_specific.OmssaParameters;
import com.compomics.util.parameters.identification.tool_specific.PNovoParameters;
import com.compomics.util.parameters.identification.tool_specific.PepnovoParameters;
import com.compomics.util.parameters.identification.tool_specific.TideParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import com.compomics.util.io.file.SerializationUtils;
import com.compomics.util.io.json.marshallers.IdentificationParametersMarshaller;
import com.compomics.util.experiment.io.parameters.DummyParameters;
import com.compomics.util.experiment.io.parameters.MarshallableParameter;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;
import com.compomics.util.parameters.identification.IdentificationParameters;
import static com.compomics.util.parameters.identification.IdentificationParameters.CURRENT_VERSION;
import com.compomics.util.parameters.identification.tool_specific.MetaMorpheusParameters;
import com.compomics.util.parameters.identification.tool_specific.SageParameters;
import com.compomics.util.parameters.identification.tool_specific.InstaNovoParameters;
import com.compomics.util.parameters.identification.tool_specific.InstaNovoPlusParameters;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * This class groups the parameters used for identification.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class SearchParameters extends ExperimentObject implements MarshallableParameter {

    /**
     * Name of the type of marshalled parameter.
     */
    private String marshallableParameterType = null;
    /**
     * Version number.
     */
    public final String version = CURRENT_VERSION;

    /**
     * Possible mass accuracy types.
     */
    public enum MassAccuracyType {

        PPM, DA;

        @Override
        public String toString() {

            switch (this) {

                case PPM:
                    return "ppm";

                case DA:
                    return "Da";

                default:
                    throw new UnsupportedOperationException(
                            "Name of tolerance type "
                            + this.name()
                            + " not implemented."
                    );
            }

        }

    };
    /**
     * Convenience array for forward ion type selection.
     */
    public static final String[] implementedForwardIons = {"a", "b", "c"};
    /**
     * Convenience array for rewind ion type selection.
     */
    public static final String[] implementedRewindIons = {"x", "y", "z"};
    /**
     * The precursor accuracy type. Default is ppm.
     */
    private MassAccuracyType precursorAccuracyType = MassAccuracyType.PPM;
    /**
     * The fragment accuracy type. Default is Da.
     */
    private MassAccuracyType fragmentAccuracyType = MassAccuracyType.PPM;
    /**
     * The precursor mass tolerance.
     */
    private double precursorTolerance = 10.0;
    /**
     * The MS2 ion tolerance.
     */
    private double fragmentIonMZTolerance = 10;
    /**
     * The expected modifications. Modified peptides will be grouped and
     * displayed according to this classification.
     */
    private ModificationParameters modificationParameters = new ModificationParameters();
    /**
     * The digestion preferences.
     */
    private DigestionParameters digestionParameters;
    /**
     * The forward ions to consider (a, b or c).
     */
    private ArrayList<Integer> forwardIons;
    /**
     * The rewind ions to consider (x, y or z).
     */
    private ArrayList<Integer> rewindIons;
    /**
     * The minimal charge searched (in absolute value).
     */
    private int minChargeSearched = 2;
    /**
     * The minimal charge searched (in absolute value).
     */
    private int maxChargeSearched = 4;
    /**
     * The minimal isotope correction.
     */
    private int minIsotopicCorrection = 0;
    /**
     * The maximal isotope correction.
     */
    private int maxIsotopicCorrection = 1;
    /**
     * Reference mass for the conversion of the fragment ion tolerance from ppm
     * to Dalton.
     */
    private double refMass = 2000.0;
    /**
     * The preferred maximal number of variable modifications.
     */
    public static final int preferredMaxVariableModifications = 6;
    /**
     * The preferred minimal number of protein sequences.
     */
    public static final int preferredMinSequences = 1000;
    /**
     * The algorithm specific parameters.
     */
    private HashMap<Integer, IdentificationAlgorithmParameter> algorithmParameters;
    /**
     * Report the flanking amino acids of a mapped peptide.
     */
    private boolean flanking = false;

    /**
     * Constructor.
     */
    public SearchParameters() {

        // Set ions to be searched by default
        forwardIons = new ArrayList<>(1);
        forwardIons.add(PeptideFragmentIon.B_ION);
        rewindIons = new ArrayList<>(1);
        rewindIons.add(PeptideFragmentIon.Y_ION);

        // Set advanced parameters
        setDefaultAdvancedSettings();

    }

    /**
     * Constructor.
     *
     * @param searchParameters the search parameter to base the search
     * parameters on.
     */
    public SearchParameters(SearchParameters searchParameters) {

        // Set values from the given parameters
        this.precursorAccuracyType = searchParameters.getPrecursorAccuracyType();
        this.fragmentAccuracyType = searchParameters.getFragmentAccuracyType();
        this.precursorTolerance = searchParameters.getPrecursorAccuracy();
        this.fragmentIonMZTolerance = searchParameters.getFragmentIonAccuracy();
        this.modificationParameters = new ModificationParameters(searchParameters.getModificationParameters());
        this.digestionParameters = searchParameters.getDigestionParameters();
        this.forwardIons = new ArrayList<>(searchParameters.getForwardIons());
        this.rewindIons = new ArrayList<>(searchParameters.getRewindIons());
        this.minChargeSearched = searchParameters.getMinChargeSearched();
        this.maxChargeSearched = searchParameters.getMaxChargeSearched();
        this.minIsotopicCorrection = searchParameters.getMinIsotopicCorrection();
        this.maxIsotopicCorrection = searchParameters.getMaxIsotopicCorrection();
        this.refMass = searchParameters.getRefMass();

        // Set advanced parameters
        setDefaultAdvancedSettings(searchParameters);

    }

    /**
     * Set the advanced settings to the default values.
     */
    public void setDefaultAdvancedSettings() {
        setDefaultAdvancedSettings(null);
    }

    /**
     * Set the advanced settings to the values in the given search parameters
     * object or to the default values if the advanced settings are not set for
     * a given advocate.
     *
     * @param searchParameters the search parameter to extract the advanced
     * settings from
     */
    public void setDefaultAdvancedSettings(SearchParameters searchParameters) {

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.omssa.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.omssa.getIndex(), new OmssaParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.omssa.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.omssa.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.xtandem.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.xtandem.getIndex(), new XtandemParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.xtandem.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.xtandem.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.msgf.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.msgf.getIndex(), new MsgfParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.msgf.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.msgf.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.msAmanda.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.msAmanda.getIndex(), new MsAmandaParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.msAmanda.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.msAmanda.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.myriMatch.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.myriMatch.getIndex(), new MyriMatchParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.myriMatch.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.myriMatch.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.comet.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.comet.getIndex(), new CometParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.comet.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.comet.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.tide.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.tide.getIndex(), new TideParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.tide.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.tide.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.andromeda.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.andromeda.getIndex(), new AndromedaParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.andromeda.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.andromeda.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.metaMorpheus.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.metaMorpheus.getIndex(), new MetaMorpheusParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.metaMorpheus.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.metaMorpheus.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.sage.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.sage.getIndex(), new SageParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.sage.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.sage.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.pepnovo.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.pepnovo.getIndex(), new PepnovoParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.pepnovo.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.pepnovo.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.direcTag.getIndex(), new DirecTagParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.direcTag.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.pNovo.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.pNovo.getIndex(), new PNovoParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.pNovo.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.pNovo.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.novor.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.novor.getIndex(), new NovorParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.novor.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.novor.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.instanovo.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.instanovo.getIndex(), new InstaNovoParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.instanovo.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.instanovo.getIndex()));
        }

        if (searchParameters == null || searchParameters.getIdentificationAlgorithmParameter(Advocate.instanovoPlus.getIndex()) == null) {
            setIdentificationAlgorithmParameter(Advocate.instanovoPlus.getIndex(), new InstaNovoPlusParameters());
        } else {
            setIdentificationAlgorithmParameter(Advocate.instanovoPlus.getIndex(), searchParameters.getIdentificationAlgorithmParameter(Advocate.instanovoPlus.getIndex()));
        }

    }

    /**
     * Returns the reference mass used to convert ppm to Da.
     *
     * @return the reference mass used to convert ppm to Da
     */
    public double getRefMass() {
        return refMass;
    }

    /**
     * Sets the reference mass used to convert ppm to Da.
     *
     * @param refMass the reference mass used to convert ppm to Da
     */
    public void setRefMass(double refMass) {
        this.refMass = refMass;
    }

    /**
     * Returns the modification settings.
     *
     * @return the modification settings
     */
    public ModificationParameters getModificationParameters() {
        return modificationParameters;
    }

    /**
     * Sets the modification settings.
     *
     * @param modificationParameters the modification settings
     */
    public void setModificationParameters(ModificationParameters modificationParameters) {
        this.modificationParameters = modificationParameters;
    }

    public void setPtmSettings(ModificationParameters modificationParameters) {
        setModificationParameters(modificationParameters);
    }

    /**
     * Returns the MS2 ion m/z tolerance.
     *
     * @return the MS2 ion m/z tolerance
     */
    public double getFragmentIonAccuracy() {
        return fragmentIonMZTolerance;
    }

    /**
     * Returns the absolute fragment ion tolerance in Dalton. If the tolerance
     * is in ppm, the internal reference mass is used.
     *
     * @return the absolute fragment ion tolerance in Dalton
     */
    public double getFragmentIonAccuracyInDaltons() {
        return getFragmentIonAccuracyInDaltons(refMass);
    }

    /**
     * Returns the absolute fragment ion tolerance in Dalton. If the tolerance
     * is in ppm, the given reference mass is used.
     *
     * @param refMass the reference mass to use for the conversion of tolerances
     * in ppm.
     *
     * @return the absolute fragment ion tolerance in Dalton
     */
    public double getFragmentIonAccuracyInDaltons(double refMass) {

        switch (fragmentAccuracyType) {

            case DA:
                return fragmentIonMZTolerance;

            case PPM:
                return fragmentIonMZTolerance * refMass / 1000000;

            default:
                throw new UnsupportedOperationException(
                        "Tolerance in "
                        + fragmentAccuracyType
                        + " not implemented."
                );

        }

    }

    /**
     * Sets the fragment ion m/z tolerance.
     *
     * @param fragmentIonMZTolerance the fragment ion m/z tolerance
     */
    public void setFragmentIonAccuracy(double fragmentIonMZTolerance) {
        this.fragmentIonMZTolerance = fragmentIonMZTolerance;
    }

    /**
     * Returns the digestion preferences.
     *
     * @return the digestion preferences
     */
    public DigestionParameters getDigestionParameters() {
        return digestionParameters;
    }

    /**
     * Sets the digestion preferences.
     *
     * @param digestionParameters the digestion preferences
     */
    public void setDigestionParameters(DigestionParameters digestionParameters) {
        this.digestionParameters = digestionParameters;
    }

    /**
     * Returns the forward ions searched as list of integers as indexed in the
     * PeptideFragmentIon class.
     *
     * @return the forward ions searched
     */
    public ArrayList<Integer> getForwardIons() {
        return forwardIons;
    }

    /**
     * Sets the forward ions searched as list of integers as indexed in the
     * PeptideFragmentIon class.
     *
     * @param forwardIons the forward ions searched
     */
    public void setForwardIons(ArrayList<Integer> forwardIons) {
        this.forwardIons = forwardIons;
    }

    /**
     * Returns the rewind ions searched as list of integers as indexed in the
     * PeptideFragmentIon class.
     *
     * @return the rewind ions searched
     */
    public ArrayList<Integer> getRewindIons() {
        return rewindIons;
    }

    /**
     * Sets the rewind ions searched as list of integers as indexed in the
     * PeptideFragmentIon class.
     *
     * @param rewindIons the rewind ions searched
     */
    public void setRewindIons(ArrayList<Integer> rewindIons) {
        this.rewindIons = rewindIons;
    }

    /**
     * Getter for the list of ion symbols used.
     *
     * @return the list of ion symbols used
     */
    public static String[] getIons() {

        String[] ions = new String[implementedForwardIons.length + implementedRewindIons.length];

        for (String forwardIon1 : implementedForwardIons) {
            ions[ions.length] = forwardIon1;
        }

        for (String rewindIon1 : implementedRewindIons) {
            ions[ions.length] = rewindIon1;
        }

        return ions;

    }

    /**
     * Returns the precursor tolerance.
     *
     * @return the precursor tolerance
     */
    public double getPrecursorAccuracy() {
        return precursorTolerance;
    }

    public double getPrecursorAccuracyDalton() {
        return precursorAccuracyType == MassAccuracyType.DA
                ? precursorTolerance
                : precursorTolerance * refMass / 1000000.0;
    }

    /**
     * Sets the precursor tolerance.
     *
     * @param precursorTolerance the precursor tolerance
     */
    public void setPrecursorAccuracy(double precursorTolerance) {
        this.precursorTolerance = precursorTolerance;
    }

    /**
     * Returns the precursor accuracy type.
     *
     * @return the precursor accuracy type
     */
    public MassAccuracyType getPrecursorAccuracyType() {
        return precursorAccuracyType;
    }

    /**
     * Sets the precursor accuracy type.
     *
     * @param precursorAccuracyType the precursor accuracy type
     */
    public void setPrecursorAccuracyType(MassAccuracyType precursorAccuracyType) {
        this.precursorAccuracyType = precursorAccuracyType;
    }

    /**
     * Returns the fragment accuracy type.
     *
     * @return the fragment accuracy type
     */
    public MassAccuracyType getFragmentAccuracyType() {
        return fragmentAccuracyType;
    }

    /**
     * Sets the fragment accuracy type.
     *
     * @param fragmentAccuracyType the fragment accuracy type
     */
    public void setFragmentAccuracyType(MassAccuracyType fragmentAccuracyType) {
        this.fragmentAccuracyType = fragmentAccuracyType;
    }

    /**
     * Returns true if the current precursor accuracy type is ppm.
     *
     * @return true if the current precursor accuracy type is ppm
     */
    public Boolean isPrecursorAccuracyTypePpm() {
        return getPrecursorAccuracyType() == MassAccuracyType.PPM;
    }

    /**
     * Returns the maximal charge searched.
     *
     * @return the maximal charge searched
     */
    public int getMaxChargeSearched() {
        return maxChargeSearched;
    }

    /**
     * Sets the maximal charge searched.
     *
     * @param maxChargeSearched the maximal charge searched
     */
    public void setMaxChargeSearched(int maxChargeSearched) {
        this.maxChargeSearched = maxChargeSearched;
    }

    /**
     * Returns the minimal charge searched.
     *
     * @return the minimal charge searched
     */
    public int getMinChargeSearched() {
        return minChargeSearched;
    }

    /**
     * Sets the minimal charge searched.
     *
     * @param minChargeSearched the minimal charge searched
     */
    public void setMinChargeSearched(int minChargeSearched) {
        this.minChargeSearched = minChargeSearched;
    }

    /**
     * Returns the algorithm specific parameters in a map: algorithm as indexed
     * in the Advocate class &gt; parameters. null if not set.
     *
     * @return the algorithm specific parameters in a map
     */
    public HashMap<Integer, IdentificationAlgorithmParameter> getAlgorithmSpecificParameters() {
        return algorithmParameters;
    }

    /**
     * Returns the algorithm specific parameters, null if not found.
     *
     * @param algorithmID the index of the search engine as indexed in the
     * Advocate class
     *
     * @return the algorithm specific parameters
     */
    public IdentificationAlgorithmParameter getIdentificationAlgorithmParameter(int algorithmID) {

        if (algorithmParameters == null) {
            return null;
        }

        return algorithmParameters.get(algorithmID);

    }

    /**
     * Adds identification algorithm specific parameters.
     *
     * @param algorithmID the algorithm id as indexed in the Advocate class
     *
     * @param identificationAlgorithmParameter the specific parameters
     */
    public void setIdentificationAlgorithmParameter(int algorithmID, IdentificationAlgorithmParameter identificationAlgorithmParameter) {

        if (algorithmParameters == null) {
            algorithmParameters = new HashMap<>();
        }

        algorithmParameters.put(algorithmID, identificationAlgorithmParameter);

    }

    /**
     * Returns the algorithms for which specific parameters are stored. Warning:
     * this does not mean that the algorithm was actually used.
     *
     * @return the algorithms for which specific parameters are stored in a set
     * of indexes as listed in the Advocate class
     */
    public Set<Integer> getAlgorithms() {

        if (algorithmParameters == null) {
            return new HashSet<>(0);
        }

        return algorithmParameters.keySet();

    }

    /**
     * Returns the minimal isotopic correction.
     *
     * @return the minimal isotopic correction
     */
    public int getMinIsotopicCorrection() {
        return minIsotopicCorrection;
    }

    /**
     * Sets the minimal isotopic correction.
     *
     * @param minIsotopicCorrection the minimal isotopic correction
     */
    public void setMinIsotopicCorrection(int minIsotopicCorrection) {
        this.minIsotopicCorrection = minIsotopicCorrection;
    }

    /**
     * Returns the maximal isotopic correction.
     *
     * @return the maximal isotopic correction
     */
    public int getMaxIsotopicCorrection() {
        return maxIsotopicCorrection;
    }

    /**
     * Sets the maximal isotopic correction.
     *
     * @param maxIsotopicCorrection the maximal isotopic correction
     */
    public void setMaxIsotopicCorrection(int maxIsotopicCorrection) {
        this.maxIsotopicCorrection = maxIsotopicCorrection;
    }

    /**
     * Loads the identification parameters from a file. If the file is an
     * identification parameters file, the search parameters are extracted.
     *
     * @param searchParametersFile the search parameter file
     *
     * @return the search parameters
     *
     * @throws IOException if an IOException occurs
     * @throws ClassNotFoundException if a ClassNotFoundException occurs
     */
    public static SearchParameters getIdentificationParameters(File searchParametersFile) throws IOException, ClassNotFoundException {

        Object savedObject;

        try {

            // Try as json file
            IdentificationParametersMarshaller jsonMarshaller = new IdentificationParametersMarshaller();
            Class expectedObjectType = DummyParameters.class;
            Object object = jsonMarshaller.fromJson(expectedObjectType, searchParametersFile);
            DummyParameters dummyParameters = (DummyParameters) object;

            if (dummyParameters.getType() == MarshallableParameter.Type.search_parameters) {
                expectedObjectType = SearchParameters.class;
                savedObject = jsonMarshaller.fromJson(expectedObjectType, searchParametersFile);
            } else if (dummyParameters.getType() == MarshallableParameter.Type.identification_parameters) {
                expectedObjectType = IdentificationParameters.class;
                savedObject = jsonMarshaller.fromJson(expectedObjectType, searchParametersFile);
            } else {
                throw new IllegalArgumentException("Parameters file " + searchParametersFile + " not recognized.");
            }

        } catch (Exception e1) {

            try {
                // try serialized java object
                savedObject = SerializationUtils.readObject(searchParametersFile);

            } catch (Exception e2) {
                e1.printStackTrace();
                e2.printStackTrace();
                throw new IllegalArgumentException("Parameters file " + searchParametersFile + " not recognized.");
            }

        }

        SearchParameters searchParameters;

        if (savedObject instanceof SearchParameters) {
            searchParameters = (SearchParameters) savedObject;
        } else if (savedObject instanceof IdentificationParameters) {
            IdentificationParameters identificationParameters = (IdentificationParameters) savedObject;
            searchParameters = identificationParameters.getSearchParameters();
        } else {
            throw new UnsupportedOperationException("Parameters of type " + savedObject.getClass() + " not supported.");
        }

        return searchParameters;

    }

    /**
     * Saves the identification parameters to a serialized file.
     *
     * @param searchParameters the identification parameters
     * @param searchParametersFile the file
     *
     * @throws IOException if an IOException occurs
     */
    public static void saveIdentificationParameters(
            SearchParameters searchParameters,
            File searchParametersFile)
            throws IOException {

        IdentificationParametersMarshaller jsonMarshaller = new IdentificationParametersMarshaller();
        searchParameters.setType();
        jsonMarshaller.saveObjectToJson(searchParameters, searchParametersFile);

    }

    /**
     * Saves the identification parameters as a human readable text file.
     *
     * @param file the file
     *
     * @throws IOException if an IOException occurs
     */
    public void saveIdentificationParametersAsTextFile(File file) throws IOException {

        try (FileWriter fw = new FileWriter(file)) {

            try (BufferedWriter bw = new BufferedWriter(fw)) {
                
                bw.write(toString());
            
            }

        }
    }

    @Override
    public String toString() {
        return toString(false);
    }

    /**
     * Returns a short description of the parameters.
     *
     * @return a short description of the parameters
     */
    public String getShortDescription() {

        SearchParameters defaultParameters = new SearchParameters();
        String newLine = System.getProperty("line.separator");
        StringBuilder output = new StringBuilder();

        if (digestionParameters != null && !DigestionParameters.getDefaultParameters().equals(digestionParameters)) {
            output.append(digestionParameters.getShortDescription());
        }

        if (modificationParameters != null) {

            ArrayList<String> modifications = modificationParameters.getFixedModifications();

            if (!modifications.isEmpty()) {
                output.append("Fixed: ");
                output.append(modifications.stream().sorted().collect(Collectors.joining(", ")));
                output.append(".").append(newLine);
            }

        }

        if (modificationParameters != null) {

            ArrayList<String> modifications = modificationParameters.getVariableModifications();

            if (!modifications.isEmpty()) {
                output.append("Variable: ");
                output.append(modifications.stream().sorted().collect(Collectors.joining(", ")));
                output.append(".").append(newLine);
            }

        }

        if (precursorTolerance != defaultParameters.getPrecursorAccuracy()
                || getPrecursorAccuracyType() != defaultParameters.getPrecursorAccuracyType()) {
            output.append("Precursor Tolerance: ").append(precursorTolerance).append(" ").append(precursorAccuracyType).append(".").append(newLine);
        }

        if (fragmentIonMZTolerance != defaultParameters.getFragmentIonAccuracy()
                || getFragmentAccuracyType() != defaultParameters.getFragmentAccuracyType()) {
            output.append("Fragment Tolerance: ").append(fragmentIonMZTolerance).append(" ").append(fragmentAccuracyType).append(".").append(newLine);
        }

        if (forwardIons != null && rewindIons != null) {

            if (!Util.sameLists(forwardIons, defaultParameters.getForwardIons())
                    || !Util.sameLists(rewindIons, defaultParameters.getRewindIons())) {

                String ions1 = forwardIons.stream()
                        .sorted()
                        .map(ion -> PeptideFragmentIon.getSubTypeAsString(ion))
                        .collect(Collectors.joining(","));

                String ions2 = rewindIons.stream()
                        .sorted()
                        .map(ion -> PeptideFragmentIon.getSubTypeAsString(ion))
                        .collect(Collectors.joining(","));

                output.append("Ion Types: ").append(ions1).append(" and ").append(ions2).append(".").append(newLine);

            }

        }

        if (minChargeSearched != defaultParameters.getMinChargeSearched()
                || maxChargeSearched != defaultParameters.getMaxChargeSearched()) {
            output.append("Charge: ").append(minChargeSearched).append("-").append(maxChargeSearched).append(".").append(newLine);
        }

        if (getMinIsotopicCorrection() != defaultParameters.getMinIsotopicCorrection()
                || getMaxIsotopicCorrection() != defaultParameters.getMaxIsotopicCorrection()) {
            output.append("Isotopic Correction: ").append(minIsotopicCorrection).append("-").append(maxIsotopicCorrection).append(".").append(newLine);
        }

        return output.toString();

    }

    /**
     * Returns the search parameters as a string.
     *
     * @param html use HTML formatting
     * @return the search parameters as a string
     */
    public String toString(boolean html) {

        String newLine;

        if (html) {
            newLine = "<br>";
        } else {
            newLine = System.getProperty("line.separator");
        }

        StringBuilder output = new StringBuilder();

        output.append("# ------------------------------------------------------------------");
        output.append(newLine);
        output.append("# General Search Parameters");
        output.append(newLine);
        output.append("# ------------------------------------------------------------------");
        output.append(newLine);
        output.append(newLine);

        if (digestionParameters.getCleavageParameter() == DigestionParameters.CleavageParameter.enzyme) {

            ArrayList<Enzyme> enzymes = digestionParameters.getEnzymes();

            for (int i = 0; i < enzymes.size(); i++) {

                Enzyme tempEnzyme = enzymes.get(i);
                String enzymeName = tempEnzyme.getName();
                output.append("ENZYME").append(i).append("=");
                output.append(enzymeName).append(", ").append(digestionParameters.getSpecificity(enzymeName));
                int nmc = digestionParameters.getnMissedCleavages(enzymeName);
                output.append(", ").append(nmc).append(" missed cleavages");
                output.append(newLine);

            }

        } else {
            output.append("ENZYME").append("=").append(digestionParameters.getCleavageParameter().name);
        }

        output.append("FIXED_MODIFICATIONS=");

        if (modificationParameters != null) {
            ArrayList<String> modifications = modificationParameters.getFixedModifications();
            output.append(modifications.stream().sorted().collect(Collectors.joining(",")));
        }

        output.append(newLine);

        output.append("VARIABLE_MODIFICATIONS=");

        if (modificationParameters != null) {
            ArrayList<String> modifications = modificationParameters.getVariableModifications();
            output.append(modifications.stream().sorted().collect(Collectors.joining(",")));
        }

        output.append(newLine);
        output.append("REFINEMENT_FIXED_MODIFICATIONS=");

        if (modificationParameters != null && modificationParameters.getRefinementFixedModifications() != null) {
            ArrayList<String> modifications = modificationParameters.getRefinementFixedModifications();
            output.append(modifications.stream().sorted().collect(Collectors.joining(",")));
        }

        output.append(newLine);
        output.append("REFINEMENT_VARIABLE_MODIFICATIONS=");

        if (modificationParameters != null && modificationParameters.getRefinementVariableModifications() != null) {
            ArrayList<String> modifications = modificationParameters.getRefinementVariableModifications();
            output.append(modifications.stream().sorted().collect(Collectors.joining(",")));
        }

        output.append(newLine);
        output.append("PRECURSOR_MASS_TOLERANCE=");
        output.append(precursorTolerance);
        output.append(newLine);
        output.append("PRECURSOR_MASS_TOLERANCE_UNIT=");

        if (getPrecursorAccuracyType() == MassAccuracyType.PPM) {
            output.append("ppm");
        } else {
            output.append("Da");
        }

        output.append(newLine);
        output.append("FRAGMENT_MASS_TOLERANCE=");
        output.append(fragmentIonMZTolerance);
        output.append(newLine);
        output.append("FRAGMENT_MASS_TOLERANCE_UNIT=");

        if (getFragmentAccuracyType() == MassAccuracyType.PPM) {
            output.append("ppm");
        } else {
            output.append("Da");
        }

        output.append(newLine);
        output.append("PPM_TO_DA_CONVERSION_REF_MASS=");
        output.append(getRefMass());
        output.append(newLine);
        output.append("FORWARD_FRAGMENT_ION_TYPE=");

        String ions1 = forwardIons.stream()
                .sorted()
                .map(ion -> PeptideFragmentIon.getSubTypeAsString(ion))
                .collect(Collectors.joining(","));

        output.append(ions1);
        output.append(newLine);
        output.append("FRAGMENT_ION_TYPE_2=");

        String ions2 = rewindIons.stream()
                .sorted()
                .map(ion -> PeptideFragmentIon.getSubTypeAsString(ion))
                .collect(Collectors.joining(","));

        output.append(ions2);
        output.append(newLine);

        output.append("PRECURSOR_CHARGE_LOWER_BOUND=");
        output.append(minChargeSearched);
        output.append(newLine);

        output.append("PRECURSOR_CHARGE_UPPER_BOUND=");
        output.append(maxChargeSearched);
        output.append(newLine);

        output.append("ISOTOPIC_CORRECTION_LOWER_BOUND=");
        output.append(getMinIsotopicCorrection());
        output.append(newLine);

        output.append("ISOTOPIC_CORRECTION_UPPER_BOUND=");
        output.append(getMaxIsotopicCorrection());
        output.append(newLine);

        for (int index : algorithmParameters.keySet()) {
            output.append(newLine);
            output.append(newLine);
            output.append(algorithmParameters.get(index).toString(html));
        }

        return output.toString();
    }

    /**
     * Returns true if the search parameter objects have identical settings.
     *
     * @param otherSearchParameters the parameters to compare to
     *
     * @return true if the search parameter objects have identical settings
     */
    public boolean equals(SearchParameters otherSearchParameters) {

        if (otherSearchParameters == null) {
            return false;
        }
        if (this.getPrecursorAccuracyType() != otherSearchParameters.getPrecursorAccuracyType()) {
            return false;
        }
        if (this.getPrecursorAccuracy() != otherSearchParameters.getPrecursorAccuracy()) {
            return false;
        }
        if (!this.getFragmentAccuracyType().equals(otherSearchParameters.getFragmentAccuracyType())) {
            return false;
        }
        if (this.getFragmentIonAccuracy() != otherSearchParameters.getFragmentIonAccuracy()) {
            return false;
        }
        if (this.getDigestionParameters() != null && otherSearchParameters.getDigestionParameters() == null
                || this.getDigestionParameters() == null && otherSearchParameters.getDigestionParameters() != null) {
            return false;
        }
        if (this.getDigestionParameters() != null && otherSearchParameters.getDigestionParameters() != null
                && !this.getDigestionParameters().isSameAs(otherSearchParameters.getDigestionParameters())) {
            return false;
        }
        if (forwardIons != null && otherSearchParameters.getForwardIons() == null
                || forwardIons == null && otherSearchParameters.getForwardIons() != null) {
            return false;
        }
        if (forwardIons != null && otherSearchParameters.getForwardIons() != null
                && !Util.sameLists(forwardIons, otherSearchParameters.getForwardIons())) {
            return false;
        }
        if (rewindIons != null && otherSearchParameters.getRewindIons() == null
                || rewindIons == null && otherSearchParameters.getRewindIons() != null) {
            return false;
        }
        if (rewindIons != null && otherSearchParameters.getRewindIons() != null
                && !Util.sameLists(rewindIons, otherSearchParameters.getRewindIons())) {
            return false;
        }
        if (this.getMinChargeSearched() != otherSearchParameters.getMinChargeSearched()) {
            return false;
        }
        if (this.getMaxChargeSearched() != otherSearchParameters.getMaxChargeSearched()) {
            return false;
        }
        if (this.getMinIsotopicCorrection() != otherSearchParameters.getMinIsotopicCorrection()) {
            return false;
        }
        if (this.getMaxIsotopicCorrection() != otherSearchParameters.getMaxIsotopicCorrection()) {
            return false;
        }
        if (!this.getModificationParameters().equals(otherSearchParameters.getModificationParameters())) {
            return false;
        }

        if (this.getAlgorithms().size() != otherSearchParameters.getAlgorithms().size()) {
            return false;
        }

        for (int se : getAlgorithms()) {

            IdentificationAlgorithmParameter otherParameter = otherSearchParameters.getIdentificationAlgorithmParameter(se);

            if (otherParameter == null) {
                return false;
            }

            IdentificationAlgorithmParameter thisParameter = getIdentificationAlgorithmParameter(se);

            if (!otherParameter.equals(thisParameter)) {
                return false;
            }
        }

        return true;

    }

    /**
     * Set whether the flanking amino acids of a mapped peptide are to be
     * reported.
     *
     * @param flanking whether the flanking amino acids of a mapped peptide are
     * to be reported
     */
    public void setFlanking(boolean flanking) {
        this.flanking = flanking;
    }

    /**
     * Return whether the flanking amino acids of a mapped peptide are to be
     * reported.
     *
     * @return true if whether the flanking amino acids of a mapped peptide are
     * to be reported
     */
    public boolean getFlanking() {
        return flanking;
    }

    @Override
    public void setType() {
        marshallableParameterType = Type.search_parameters.name();
    }

    @Override
    public Type getType() {

        if (marshallableParameterType == null) {
            return null;
        }

        return Type.valueOf(marshallableParameterType);

    }
}
