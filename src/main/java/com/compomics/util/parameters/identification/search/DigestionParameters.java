package com.compomics.util.parameters.identification.search;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.aminoacids.AminoAcid;
import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.personalization.ExperimentObject;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * This class groups the parameters for the digestion of proteins.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class DigestionParameters extends ExperimentObject {

    /**
     * Enum for the different types of enzyme specificity.
     */
    public enum Specificity {

        /**
         * Specific at both termini.
         */
        specific(0, "Specific"),
        /**
         * Specific at only one of the termini.
         */
        semiSpecific(1, "Semi-Specific"),
        /**
         * Specific at the N-terminus only.
         */
        specificNTermOnly(2, "N-term Specific"),
        /**
         * Specific at the C-terminus only.
         */
        specificCTermOnly(3, "C-term Specific");

        /**
         * The index.
         */
        public final int index;
        /**
         * The name.
         */
        public final String name;

        /**
         * Constructor.
         *
         * @param index the index as integer
         * @param name the name
         */
        private Specificity(int index, String name) {
            this.index = index;
            this.name = name;
        }

        /**
         * Returns the specificity of the given index.
         *
         * @param index the index of the specificity
         *
         * @return the corresponding specificity
         */
        public static Specificity getSpecificity(int index) {

            for (Specificity specificity : values()) {

                if (specificity.index == index) {
                    return specificity;
                }

            }

            throw new IllegalArgumentException("No specificity found for index " + index + ".");
        }

        /**
         * Returns the different options as command line description.
         *
         * @return the different options as command line description
         */
        public static String getCommandLineDescription() {

            StringBuilder stringBuilder = new StringBuilder();

            for (Specificity specificity : values()) {

                if (stringBuilder.length() > 0) {
                    stringBuilder.append(", ");
                }

                stringBuilder.append(specificity.index).append(": ").append(specificity.name);

            }

            return stringBuilder.toString();
        }

        @Override
        public String toString() {
            return name;
        }
    }

    public static DigestionParameters getDefaultPreferences() {
        return getDefaultParameters();
    }

    /**
     * Enum for the different types of digestion.
     */
    public enum CleavageParameter {

        /**
         * Digestion with an enzyme.
         */
        enzyme(0, "Enzyme"),
        /**
         * Unspecific digestion.
         */
        unSpecific(1, "Unspecific"),
        /**
         * Whole protein, no digestion.
         */
        wholeProtein(2, "Whole Protein");

        /**
         * The index.
         */
        public final int index;
        /**
         * The name.
         */
        public final String name;

        /**
         * Constructor.
         *
         * @param index the index as integer
         * @param name the name
         */
        private CleavageParameter(int index, String name) {
            this.index = index;
            this.name = name;
        }

        /**
         * Returns the cleavage parameter of the given index.
         *
         * @param index the index of the cleavage parameter
         *
         * @return the corresponding cleavage parameter
         */
        public static CleavageParameter getCleavageParameters(int index) {

            for (CleavageParameter cleavageParameter : values()) {

                if (cleavageParameter.index == index) {
                    return cleavageParameter;
                }

            }

            throw new IllegalArgumentException("No cleavage parameter found for index " + index + ".");
        }

        /**
         * Returns the different options as command line description.
         *
         * @return the different options as command line description
         */
        public static String getCommandLineDescription() {

            StringBuilder stringBuilder = new StringBuilder();

            for (CleavageParameter cleavageParameter : values()) {

                if (stringBuilder.length() > 0) {
                    stringBuilder.append(", ");
                }

                stringBuilder.append(cleavageParameter.index).append(": ").append(cleavageParameter.name);

            }

            return stringBuilder.toString();
        }

        @Override
        public String toString() {
            return name;
        }
    }

    public enum CleavagePreference {
        enzyme,
        unSpecific,
        wholeProtein
    }

    /**
     * Boolean indicating whether the sample was not digested.
     */
    private CleavageParameter cleavageParameter;
    /**
     * List of enzyme used.
     */
    private ArrayList<Enzyme> enzymes;
    /**
     * Number of allowed missed cleavages.
     */
    private HashMap<String, Integer> nMissedCleavages;
    /**
     * The specificity of the enzyme.
     */
    private HashMap<String, Specificity> specificity;

    /**
     * Constructor for empty parameters.
     */
    public DigestionParameters() {
    }

    /**
     * Clones the given parameters.
     *
     * @param digestionParameters the parameters to clone
     *
     * @return a new object containing the same parameters
     */
    public static DigestionParameters clone(DigestionParameters digestionParameters) {

        DigestionParameters clone = new DigestionParameters();
        clone.setCleavageParameter(digestionParameters.getCleavageParameter());

        if (digestionParameters.getCleavageParameter() == DigestionParameters.CleavageParameter.enzyme) {

            for (Enzyme enzyme : digestionParameters.getEnzymes()) {

                clone.addEnzyme(enzyme);
                String enzymeName = enzyme.getName();
                clone.setSpecificity(enzymeName, digestionParameters.getSpecificity(enzymeName));
                clone.setnMissedCleavages(enzymeName, digestionParameters.getnMissedCleavages(enzymeName));

            }

        }

        return clone;
    }

    /**
     * Returns default digestion parameters. Trypsin specific with two missed
     * cleavages.
     *
     * @return default digestion parameters
     */
    public static DigestionParameters getDefaultParameters() {

        DigestionParameters digestionParameters = new DigestionParameters();
        digestionParameters.setCleavageParameter(CleavageParameter.enzyme);
        String enzymeName = "Trypsin";
        Enzyme trypsin = EnzymeFactory.getInstance().getEnzyme(enzymeName);

        if (trypsin != null) { // required if the user deleted trypsin from the factory
            digestionParameters.addEnzyme(trypsin);
            digestionParameters.setnMissedCleavages(enzymeName, 2);
        }

        return digestionParameters;
    }

    /**
     * Returns a boolean indicating whether enzyme settings were set.
     *
     * @return a boolean indicating whether enzyme settings were set
     */
    public boolean hasEnzymes() {

        return enzymes != null && !enzymes.isEmpty();

    }

    /**
     * Returns the enzymes used for digestion in a list.
     *
     * @return the enzymes used for digestion in a list
     */
    public ArrayList<Enzyme> getEnzymes() {

        return enzymes;

    }

    /**
     * Sets the enzymes used for digestion.
     *
     * @param enzymes the enzymes used for digestion in a list
     */
    public void setEnzymes(ArrayList<Enzyme> enzymes) {

        this.enzymes = enzymes;

    }

    /**
     * Adds an enzyme. The specificity of the enzyme is set by default to
     * specific and the number of allowed missed cleavages to 0.
     *
     * @param enzyme an enzyme used for digestion.
     */
    public void addEnzyme(Enzyme enzyme) {

        if (enzymes == null) {
            enzymes = new ArrayList<>(1);
        }

        enzymes.add(enzyme);
        setSpecificity(enzyme.getName(), Specificity.specific);
        setnMissedCleavages(enzyme.getName(), 0);

    }

    /**
     * Clears the parameters.
     */
    public void clear() {

        cleavageParameter = null;
        enzymes = null;
        nMissedCleavages = null;
        specificity = null;

    }

    /**
     * Clears the enzymes set including specificity and missed cleavages.
     */
    public void clearEnzymes() {

        enzymes = null;
        nMissedCleavages = null;
        specificity = null;

    }

    /**
     * Returns the number of allowed missed cleavages for the given enzyme. Null
     * if not set.
     *
     * @param enzymeName the name of the enzyme
     *
     * @return the number of allowed missed cleavages
     */
    public Integer getnMissedCleavages(String enzymeName) {

        if (nMissedCleavages == null) {
            return null;
        }

        return nMissedCleavages.get(enzymeName);

    }

    /**
     * Sets the number of allowed missed cleavages.
     *
     * @param enzymeName the name of the enzyme
     * @param enzymeMissedCleavages the number of allowed missed cleavages
     */
    public void setnMissedCleavages(String enzymeName, int enzymeMissedCleavages) {

        if (nMissedCleavages == null) {
            nMissedCleavages = new HashMap<>(1);
        }

        nMissedCleavages.put(enzymeName, enzymeMissedCleavages);

    }

    /**
     * Returns the expected specificity of the given enzyme. Null if not set.
     *
     * @param enzymeName the name of the enzyme
     *
     * @return the specificity
     */
    public Specificity getSpecificity(String enzymeName) {

        if (specificity == null) {
            return null;
        }

        return specificity.get(enzymeName);

    }

    /**
     * Sets the expected specificity of the enzyme.
     *
     * @param enzymeName the name of the enzyme
     * @param enzymeSpecificity the expected specificity of the enzyme
     */
    public void setSpecificity(String enzymeName, Specificity enzymeSpecificity) {

        if (specificity == null) {
            specificity = new HashMap<>(1);
        }

        specificity.put(enzymeName, enzymeSpecificity);

    }

    /**
     * Returns the cleavage parameters.
     *
     * @return the cleavage parameters
     */
    public CleavageParameter getCleavageParameter() {

        return cleavageParameter;

    }

    public CleavagePreference getCleavagePreference() {
        return CleavagePreference.valueOf(cleavageParameter.name());
    }

    /**
     * Sets the cleavage parameters.
     *
     * @param cleavageParameter the cleavage parameters
     */
    public void setCleavageParameter(CleavageParameter cleavageParameter) {

        this.cleavageParameter = cleavageParameter;

    }

    /**
     * Returns a short description of the parameters.
     *
     * @return a short description of the parameters
     */
    public String getShortDescription() {

        DigestionParameters defaultParameters = DigestionParameters.getDefaultParameters();
        StringBuilder stringBuilder = new StringBuilder();
        String newLine = System.getProperty("line.separator");

        if (!defaultParameters.isSameAs(this)) {

            stringBuilder.append("Digestion: ");

            switch (cleavageParameter) { // @TODO: can be null..?

                case wholeProtein:
                    stringBuilder.append("Whole Protein").append(newLine);
                    break;

                case unSpecific:
                    stringBuilder.append("Unspecific").append(newLine);
                    break;

                case enzyme:
                    for (Enzyme enzyme : enzymes) {

                        if (stringBuilder.length() > 0) {
                            stringBuilder.append(newLine);
                        }

                        String enzymeName = enzyme.getName();
                        stringBuilder.append(enzymeName).append(", ");
                    }
                    
                    // @TODO: support enzyme-specific specificity and missed cleavages?
                    Enzyme firstEnzyme = enzymes.get(0);
                    
                    stringBuilder.append(getSpecificity(firstEnzyme.getName()));
                    
                    Integer nmc = getnMissedCleavages(firstEnzyme.getName());

                    if (nmc != null) {
                        stringBuilder.append(", ").append(nmc).append(" missed cleavages");
                    }
                    
                    break;

                default:
                    throw new UnsupportedOperationException(
                            "Description not implemented for cleavage parameter "
                            + cleavageParameter
                            + "."
                    );
            }

            stringBuilder.append(".").append(newLine);
        }

        return stringBuilder.toString();
    }

    /**
     * Returns a boolean indicating whether these digestion parameters are the
     * same as the given other parameters.
     *
     * @param otherDigestionParameters the other digestion parameters
     *
     * @return a boolean indicating whether these digestion parameters are the
     * same as the given other parameters
     */
    public boolean isSameAs(DigestionParameters otherDigestionParameters) {

        if (cleavageParameter != otherDigestionParameters.getCleavageParameter()) {
            return false;
        }

        ArrayList<Enzyme> otherEnzymes = otherDigestionParameters.getEnzymes();

        if ((enzymes != null && otherEnzymes == null)
                || (enzymes == null && otherEnzymes != null)) {
            return false;
        }

        if (enzymes != null && otherEnzymes != null) {

            if (enzymes.size() != otherEnzymes.size()) {
                return false;
            }

            ArrayList<String> enzymeNames = new ArrayList<>(enzymes.size());

            for (Enzyme enzyme : enzymes) {
                enzymeNames.add(enzyme.getName());
            }

            ArrayList<String> otherNames = new ArrayList<>(otherEnzymes.size());

            for (Enzyme enzyme : otherEnzymes) {
                otherNames.add(enzyme.getName());
            }

            if (!Util.sameLists(enzymeNames, otherNames)) {
                return false;
            }

            for (String enzymeName : enzymeNames) {

                if (getSpecificity(enzymeName) != otherDigestionParameters.getSpecificity(enzymeName)
                        || !getnMissedCleavages(enzymeName).equals(otherDigestionParameters.getnMissedCleavages(enzymeName))) {
                    return false;
                }

            }
        }

        return true;
    }

    /**
     * Get the X! Tandem enzyme format.
     *
     * @return the enzyme X! Tandem format as String
     */
    public String getXTandemFormat() {

        switch (cleavageParameter) {

            case wholeProtein:
                return "";

            case unSpecific:
                return "[X]|[X]";

            case enzyme:
                StringBuilder result = new StringBuilder();

                for (Enzyme enzyme : enzymes) {

                    if (result.length() > 0) {
                        result.append(",");
                    }

                    if (!enzyme.getAminoAcidBefore().isEmpty()) {

                        result.append("[");

                        for (Character aa : enzyme.getAminoAcidBefore()) {
                            result.append(aa);
                        }

                        result.append("]");
                    }

                    if (!enzyme.getRestrictionBefore().isEmpty()) {

                        result.append("{");

                        for (Character aa : enzyme.getRestrictionBefore()) {
                            result.append(aa);
                        }

                        result.append("}");
                    }

                    if (enzyme.getAminoAcidBefore().isEmpty()
                            && enzyme.getRestrictionBefore().isEmpty()) {
                        result.append("[X]");
                    }

                    result.append("|");

                    if (!enzyme.getAminoAcidAfter().isEmpty()) {

                        result.append("[");

                        for (Character aa : enzyme.getAminoAcidAfter()) {
                            result.append(aa);
                        }

                        result.append("]");
                    }

                    if (!enzyme.getRestrictionAfter().isEmpty()) {

                        result.append("{");

                        for (Character aa : enzyme.getRestrictionAfter()) {
                            result.append(aa);
                        }

                        result.append("}");
                    }

                    if (enzyme.getAminoAcidAfter().isEmpty()
                            && enzyme.getRestrictionAfter().isEmpty()) {
                        result.append("[X]");
                    }
                }
                
                return result.toString();

            default:
                throw new UnsupportedOperationException(
                        "X! Tandem format not implemented for cleavage parameter "
                        + cleavageParameter
                        + "."
                );
        }
    }

    /**
     * Get the MyriMatch enzyme format. In case multiple enzymes are present all
     * possible cleavage sites will be included.
     *
     * @return the enzyme MyriMatch format as String
     */
    public String getMyriMatchFormat() {

        // example: trypsin corresponds to "[|R|K . . ]"
        // details: http://htmlpreview.github.io/?https://github.com/ProteoWizard/pwiz/blob/master/pwiz_tools/Bumbershoot/myrimatch/doc/index.html
        String result = "[";

        for (Enzyme enzyme : enzymes) {

            if (enzyme.getAminoAcidBefore().size() > 0) {

                for (Character aa : enzyme.getAminoAcidBefore()) {
                    result += "|" + aa;
                }

                result += " ";

            } else {

                result += " ";

            }
        }

        HashSet<Character> commonRestrictionAfter;

        if (enzymes.size() == 1) {

            commonRestrictionAfter = enzymes.get(0).getRestrictionAfter();

        } else {

            commonRestrictionAfter = new HashSet<>();

            for (Character aa : enzymes.get(0).getRestrictionAfter()) {

                boolean missing = false;

                for (Enzyme enzyme : enzymes) {

                    if (!enzyme.getRestrictionAfter().contains(aa)) {
                        missing = true;
                        break;
                    }
                }

                if (!missing) {
                    commonRestrictionAfter.add(aa);
                }
            }
        }

        if (commonRestrictionAfter.size() > 0) {

            String temp = "";

            for (Character aa : AminoAcid.getUniqueAminoAcids()) {

                if (!commonRestrictionAfter.contains(aa)) {

                    if (!temp.isEmpty()) {
                        temp += "|";
                    }

                    temp += aa;
                }
            }

            result += temp + " ";

        } else {
            result += ". ";
        }

        HashSet<Character> commonRestrictionBefore;

        if (enzymes.size() == 1) {

            commonRestrictionBefore = enzymes.get(0).getRestrictionBefore();

        } else {

            commonRestrictionBefore = new HashSet<>();

            for (Character aa : enzymes.get(0).getRestrictionBefore()) {

                boolean missing = false;

                for (Enzyme enzyme : enzymes) {

                    if (!enzyme.getRestrictionBefore().contains(aa)) {
                        missing = true;
                        break;
                    }

                }

                if (!missing) {
                    commonRestrictionBefore.add(aa);
                }
            }
        }

        if (commonRestrictionBefore.size() > 0) {

            String temp = "";

            for (Character aa : AminoAcid.getUniqueAminoAcids()) {

                if (!commonRestrictionBefore.contains(aa)) {

                    if (!temp.isEmpty()) {
                        temp += "|";
                    }

                    temp += aa;
                }
            }

            result += temp + " ";

        } else {
            result += ". ";
        }

        for (Enzyme enzyme : enzymes) {

            if (enzyme.getAminoAcidAfter().size() > 0) {

                String temp = "";

                for (Character aa : enzyme.getAminoAcidAfter()) {

                    if (!temp.isEmpty()) {

                        temp += "|";
                    }

                    temp += aa;
                }

                result += temp + "|";
            }
        }

        return result + "]";
    }
}
