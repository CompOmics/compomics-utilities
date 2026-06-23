package com.compomics.util.experiment.biology.modifications;

import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidPattern;
import com.compomics.util.experiment.biology.atoms.Atom;
import com.compomics.util.experiment.biology.atoms.AtomChain;
import com.compomics.util.experiment.biology.atoms.AtomImpl;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import com.compomics.util.experiment.biology.ions.impl.ReporterIon;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.io.json.JsonMarshaller;
import com.compomics.util.parameters.identification.tool_specific.MetaMorpheusParameters;
import com.compomics.util.pride.CvTerm;
import java.awt.Color;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.stream.Collectors;

/**
 * This factory will load Modification from an XML file and provide them on
 * demand as standard class.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class ModificationFactory implements ModificationProvider {
    
    /**
     * Instance of the factory.
     */
    private static ModificationFactory instance = null;
    /**
     * The folder containing the Modification factory.
     */
    private static String SERIALIZATION_FILE_FOLDER = System.getProperty("user.home") + "/.compomics";
    /**
     * The name of the Modification factory back-up file. The version number
     * follows the one of utilities.
     */
    private static final String SERIALIZATION_FILE_NAME = "modificationFactory-5.0.42.json";
    /**
     * A map linking indexes with modifications.
     */
    private final HashMap<String, Modification> modificationMap = new HashMap<>();
    /**
     * List of the indexes of default modifications.
     */
    private final ArrayList<String> defaultMods = new ArrayList<>();
    /**
     * List of the indexes of user modifications.
     */
    private final ArrayList<String> userMods = new ArrayList<>();
    /**
     * Mapping of the expected modification names to the color used.
     */
    private final HashMap<String, Integer> userColors = new HashMap<>();
    /**
     * Map of modification names mapping to a given PSI-MOD accession number
     * (key provided without the "MOD:" part).
     */
    private final HashMap<String, ArrayList<String>> psiModMap = new HashMap<>();
    /**
     * Suffix for the modification clone targeting a single amino acid instead
     * of a pattern.
     */
    public static final String SINGLE_AA_SUFFIX = "|single_aa";
    /**
     * Set to true if the default mods are sorted alphabetically.
     */
    public boolean defaultModsSorted = false;
    /**
     * Set to true if the users mods are sorted alphabetically.
     */
    public boolean usersModsSorted = false;

    /**
     * Constructor for the factory.
     */
    private ModificationFactory() {
        defaultModsSorted = false;
        setDefaultModifications();
    }

    /**
     * Static method to get the instance of the factory.
     *
     * @return the instance of the factory
     */
    public static ModificationFactory getInstance() {
        if (instance == null) {
            try {
                File savedFile = new File(SERIALIZATION_FILE_FOLDER, SERIALIZATION_FILE_NAME);
                instance = loadFromFile(savedFile);
                instance.checkUserModifications();
            } catch (Exception e) {
                instance = new ModificationFactory();
            }
        }
        return instance;
    }

    /**
     * Loads an enzyme factory from a file. The file must be an export of the
     * factory in the json format.
     *
     * @param file the file to load
     *
     * @return the enzyme factory saved in file
     *
     * @throws IOException exception thrown whenever an error occurred while
     * loading the file
     */
    public static ModificationFactory loadFromFile(File file) throws IOException {
        JsonMarshaller jsonMarshaller = new JsonMarshaller();
        ModificationFactory result = (ModificationFactory) jsonMarshaller.fromJson(ModificationFactory.class, file);
        return result;
    }

    /**
     * Saves a Modification factory to a file.
     *
     * @param modificationFactory the Modification factory to save
     * @param file the file where to save
     *
     * @throws IOException exception thrown whenever an error occurred while
     * saving the file
     */
    public static void saveToFile(ModificationFactory modificationFactory, File file) throws IOException {
        JsonMarshaller jsonMarshaller = new JsonMarshaller();
        jsonMarshaller.saveObjectToJson(modificationFactory, file);
    }

    /**
     * Add neutral losses and reporter ions for the user Modifications.
     */
    private void checkUserModifications() {
        for (String tempUserMod : getUserModifications()) {

            Modification modification = getModification(tempUserMod);

            if (!modification.getNeutralLosses().isEmpty()) {
                for (NeutralLoss neutralLoss : modification.getNeutralLosses()) {
                    if (NeutralLoss.getNeutralLoss(neutralLoss.name) == null) {
                        NeutralLoss.addNeutralLoss(neutralLoss);
                    }
                }
            }
            if (!modification.getReporterIons().isEmpty()) {
                for (ReporterIon reporterIon : modification.getReporterIons()) {
                    if (ReporterIon.getReporterIon(reporterIon.getName()) == null) {
                        ReporterIon.addReporterIon(reporterIon);
                    }
                }
            }
        }
    }

    /**
     * Clears the factory getInstance() needs to be called afterwards.
     */
    public void clearFactory() {
        instance = new ModificationFactory();
    }

    /**
     * Reloads the factory getInstance() needs to be called afterwards.
     */
    public void reloadFactory() {
        instance = null;
    }

    /**
     * Saves the factory in the user folder.
     *
     * @throws IOException exception thrown whenever an error occurred while
     * saving the modificationFactory
     */
    public void saveFactory() throws IOException {
        File factoryFile = new File(SERIALIZATION_FILE_FOLDER, SERIALIZATION_FILE_NAME);
        if (!factoryFile.getParentFile().exists()) {
            factoryFile.getParentFile().mkdir();
        }
        saveToFile(instance, factoryFile);
    }

    /**
     * Returns a clone of the given Modification targeting a single amino acid
     * instead of a pattern.
     *
     * @param modification the modification of interest
     *
     * @return a clone of the given Modification targeting a single amino acid
     * instead of a pattern
     */
    public static Modification getSingleAAModification(Modification modification) {
        if (!modification.isStandardSearch()) {
            return new Modification(modification.getModificationType(), modification.getShortName(),
                    modification.getName() + SINGLE_AA_SUFFIX, modification.getAtomChainAdded(),
                    modification.getAtomChainRemoved(), modification.getPattern().getStandardSearchPattern(),
                    modification.getCategory());
        } else {
            return modification;
        }
    }

    /**
     * Returns a clone of the given Modification targeting a single amino acid
     * instead of a pattern.
     *
     * @param modificationName the name of the modification of interest
     *
     * @return a clone of the given Modification targeting a single amino acid
     * instead of a pattern
     */
    public Modification getSingleAAModification(String modificationName) {
        Modification modification = getModification(modificationName);
        return getSingleAAModification(modification);
    }

    /**
     * Adds a new user modification.
     *
     * @param modification the new modification to add
     */
    public void addUserModification(Modification modification) {

        String modn_proteiname = modification.getName();
        modificationMap.put(modn_proteiname, modification);
        if (!userMods.contains(modn_proteiname)) {
            userMods.add(modn_proteiname);
        } else {
            userMods.set(userMods.indexOf(modn_proteiname), modn_proteiname);
        }
        usersModsSorted = false;

        // add the neutral losses and reporter ions
        for (NeutralLoss neutralLoss : modification.getNeutralLosses()) {
            if (NeutralLoss.getNeutralLoss(neutralLoss.name) == null) {
                NeutralLoss.addNeutralLoss(neutralLoss);
            }
        }
        for (ReporterIon reporterIon : modification.getReporterIons()) {
            if (ReporterIon.getReporterIon(reporterIon.getName()) == null) {
                ReporterIon.addReporterIon(reporterIon);
            }
        }
    }

    /**
     * Removes a user Modification.
     *
     * @param modificationName the name of the Modification to remove
     */
    public void removeUserPtm(String modificationName) {
        if (defaultMods.contains(modificationName)) {
            throw new IllegalArgumentException("Impossible to remove default modification " + modificationName);
        }
        modificationMap.remove(modificationName);
        userMods.remove(modificationName);
    }

    @Override
    public Modification getModification(String name) {
        return modificationMap.get(name);
    }

    /**
     * Returns a boolean indicating whether the Modification is loaded in the
     * factory.
     *
     * @param name the name of the Modification
     * @return a boolean indicating whether the Modification is loaded in the
     * factory
     */
    public boolean containsModification(String name) {
        return modificationMap.containsKey(name);
    }

    public boolean containsPTM(String name) {
        return containsModification(name);
    }

    /**
     * Returns the names of the default modifications.
     *
     * @return the names of the default modifications
     */
    public ArrayList<String> getDefaultModifications() {
        return defaultMods;
    }

    /**
     * Returns the alphabetically (case insensitive) ordered names of the
     * default modifications.
     *
     * @return the alphabetically ordered names of the default modifications
     */
    public ArrayList<String> getDefaultModificationsOrdered() {
        if (!defaultModsSorted) {
            Collections.sort(defaultMods, String.CASE_INSENSITIVE_ORDER);
            defaultModsSorted = true;
        }
        return defaultMods;
    }

    /**
     * Returns the names of the user defined modifications.
     *
     * @return the names of the user defined modifications
     */
    public ArrayList<String> getUserModifications() {
        return userMods;
    }

    /**
     * Returns the alphabetically (case insensitive) ordered names of the user
     * defined modifications.
     *
     * @return the alphabetically ordered names of the user defined
     * modifications
     */
    public ArrayList<String> getUserModificationsOrdered() {
        if (!usersModsSorted) {
            Collections.sort(userMods, String.CASE_INSENSITIVE_ORDER);
            usersModsSorted = true;
        }
        return userMods;
    }

    /**
     * Returns the names of all imported Modifications.
     *
     * @return the names of all imported Modifications
     */
    public ArrayList<String> getModifications() {
        return new ArrayList<>(modificationMap.keySet());
    }

    public ArrayList<String> getPTMs() {
        return getModifications();
    }

    /**
     * Convenience method returning a boolean indicating whether a Modification
     * is user defined or default.
     *
     * @param modificationName the name of the Modification
     * @return boolean indicating whether a Modification is user defined
     */
    public boolean isUserDefined(String modificationName) {
        return !defaultMods.contains(modificationName);
    }

    /**
     * Verifies that the modifications backed-up in the search parameters are
     * loaded and alerts the user in case conflicts are found.
     *
     * @param searchParameters the search parameters to load
     * @param overwrite if true, overwrite the modification
     * @return returns a list of modifications already loaded which should be
     * checked.
     */
    public ArrayList<String> loadBackedUpModifications(SearchParameters searchParameters, boolean overwrite) {
        ModificationParameters modificationProfile = searchParameters.getModificationParameters();
        ArrayList<String> toCheck = new ArrayList<>();
        for (String newModificationName : modificationProfile.getBackedUpModifications().keySet()) {
            if (containsModification(newModificationName)) {
                Modification oldModification = getModification(newModificationName);
                Modification newModification = modificationProfile.getModification(newModificationName);
                if (!oldModification.isSameAs(newModification)) {
                    toCheck.add(newModificationName);
                    if (overwrite) {
                        modificationMap.put(newModificationName, newModification);
                        for (NeutralLoss neutralLoss : newModification.getNeutralLosses()) {
                            NeutralLoss.addNeutralLoss(neutralLoss);
                        }
                    }
                }
            } else {
                Modification modification = modificationProfile.getModification(newModificationName);
                addUserModification(modification);
                for (NeutralLoss neutralLoss : modification.getNeutralLosses()) {
                    NeutralLoss.addNeutralLoss(neutralLoss);
                }
            }
        }
        return toCheck;
    }

    /**
     * Returns the color used to code the given modification.
     *
     * @param modification the name of the given expected modification
     * @return the corresponding color
     */
    public int getColor(String modification) {
        if (!userColors.containsKey(modification)) {
            setColor(modification, getDefaultColor(modification));
        }
        return userColors.get(modification);
    }

    /**
     * Sets a new color for the given expected modification.
     *
     * @param expectedModification the name of the expected modification
     * @param color the new color
     */
    public void setColor(String expectedModification, int color) {
        userColors.put(expectedModification, color);
    }

    public void setColor(String expectedModification, Color color) {
        setColor(expectedModification, color.getRGB());
    }

    public void addFixedModifications(
            ModificationParameters modificationParameters,
            com.compomics.util.experiment.identification.amino_acid_tags.Tag tag,
            com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters sequenceMatchingParameters
    ) {
        // Fixed modifications are resolved on demand in Utilities 5.
    }

    /**
     * Returns a default color based on the modification name.
     *
     * @param modification the name of the modification
     * @return a default color.
     */
    public static int getDefaultColor(String modification) {
        if (modification.contains("no modification")) {
            return Color.LIGHT_GRAY.getRGB();
        } else if (modification.toLowerCase().contains("phospho")) {
            return Color.RED.getRGB();
        } else if (modification.toLowerCase().contains("pyro")) {
            return new Color(255, 102, 51).getRGB();
        } else if (modification.toLowerCase().contains("ox")) {
            return Color.BLUE.getRGB();
        } else if (modification.toLowerCase().contains("itraq")) {
            return Color.ORANGE.getRGB();
        } else if (modification.toLowerCase().contains("tmt")) {
            return Color.ORANGE.getRGB();
        } else if (modification.toLowerCase().contains("carbamido")) {
            return Color.LIGHT_GRAY.getRGB();
        } else if (modification.toLowerCase().contains("ace")) {
            return new Color(153, 153, 0).getRGB();
        } else if (modification.toLowerCase().contains("glyco")) {
            return Color.MAGENTA.getRGB();
        } else {
            float r = (float) Math.random();
            float g = (float) Math.random();
            float b = (float) Math.random();
            return new Color(r, g, b).getRGB();
        }
    }

    /**
     * Tries to convert a PRIDE Modification to utilities Modification name, and
     * add it to the modification profile. Unknown Modifications are added to
     * the unknown Modifications arraylist.
     *
     * @param pridePtmName the PRIDE Modification name
     * @param modProfile the modification profile to add the Modifications to
     * @param unknownPtms the list of unknown ModificationS, updated during this
     * method
     * @param isFixed if true, the Modification will be added as a fixed
     * modification
     * @return a pride parameters report as a string (for use in PRIDE Reshake)
     */
    public String convertPridePtm(String pridePtmName, ModificationParameters modProfile, ArrayList<String> unknownPtms, boolean isFixed) {

        String prideParametersReport = "";

        // special cases for when multiple modifications are needed
        if (pridePtmName.equalsIgnoreCase("iTRAQ4plex")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex reporter+balance reagent N-acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex-114 reporter+balance reagent N6-acylated lysine")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex-114 reporter+balance reagent O4&apos;-acylated tyrosine")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex-114 reporter+balance reagent acylated N-terminal")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex-114 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex-116 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ4plex-117 reporter+balance reagent N6-acylated lysine")) {

            if (!modProfile.contains("iTRAQ 4-plex of K")) {
                modProfile.addFixedModification(getModification("iTRAQ 4-plex of K"));
                prideParametersReport += "<br>" + "iTRAQ 4-plex of K" + " (assumed fixed)";
            }
            if (!modProfile.contains("iTRAQ 4-plex of N-term")) {
                modProfile.addFixedModification(getModification("iTRAQ 4-plex of N-term"));
                prideParametersReport += "<br>" + "iTRAQ 4-plex of N-term" + " (assumed fixed)";
            }
            if (!modProfile.contains("iTRAQ 4-plex of Y")) {
                modProfile.addVariableModification(getModification("iTRAQ 4-plex of Y"));
                prideParametersReport += "<br>" + "iTRAQ 4-plex of Y" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("iTRAQ8plex")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex-113 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex-114 reporter+balance reagent N6-acylated lysine")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex-114 reporter+balance reagent O4&apos;-acylated tyrosine")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex-114 reporter+balance reagent acylated N-terminal")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex-115 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex-116 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("iTRAQ8plex:13C(6)15N(2)")) {

            if (!modProfile.contains("iTRAQ 8-plex of K")) {
                modProfile.addFixedModification(getModification("iTRAQ 8-plex of K"));
                prideParametersReport += "<br>" + "iTRAQ 8-plex of K" + " (assumed fixed)";
            }
            if (!modProfile.contains("iTRAQ 8-plex of N-term")) {
                modProfile.addFixedModification(getModification("iTRAQ 8-plex of N-term"));
                prideParametersReport += "<br>" + "iTRAQ 8-plex of N-term" + " (assumed fixed)";
            }
            if (!modProfile.contains("iTRAQ 8-plex of Y")) {
                modProfile.addVariableModification(getModification("iTRAQ 8-plex of Y"));
                prideParametersReport += "<br>" + "iTRAQ 8-plex of Y" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("TMT2plex") || pridePtmName.equalsIgnoreCase("TMTduplex")) {

            if (!modProfile.contains("TMT 2-plex of K")) {
                modProfile.addFixedModification(getModification("TMT 2-plex of K"));
                prideParametersReport += "<br>" + "TMT 2-plex of K" + " (assumed fixed)";
            }
            if (!modProfile.contains("TMT 2-plex of N-term")) {
                modProfile.addFixedModification(getModification("TMT 2-plex of N-term"));
                prideParametersReport += "<br>" + "TMT 2-plex of N-term" + " (assumed fixed)";
            }

        } else if (pridePtmName.equalsIgnoreCase("TMT6plex")
                || pridePtmName.equalsIgnoreCase("TMT6plex-126 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("TMT6plex-131 reporter+balance reagent N6-acylated lysine")) {

            if (!modProfile.contains("TMT 6-plex of K")) {
                modProfile.addFixedModification(getModification("TMT 6-plex of K"));
                prideParametersReport += "<br>" + "TMT 6-plex of K" + " (assumed fixed)";
            }
            if (!modProfile.contains("TMT 6-plex of N-term")) {
                modProfile.addFixedModification(getModification("TMT 6-plex of N-term"));
                prideParametersReport += "<br>" + "TMT 6-plex of N-term" + " (assumed fixed)";
            }

        } else if (pridePtmName.equalsIgnoreCase("TMT10plex")
                || pridePtmName.equalsIgnoreCase("TMT10plex-126 reporter+balance reagent acylated residue")
                || pridePtmName.equalsIgnoreCase("TMT10plex-131 reporter+balance reagent N6-acylated lysine")) {

            if (!modProfile.contains("TMT 10-plex of K")) {
                modProfile.addFixedModification(getModification("TMT 10-plex of K"));
                prideParametersReport += "<br>" + "TMT 10-plex of K" + " (assumed fixed)";
            }
            if (!modProfile.contains("TMT 10-plex of N-term")) {
                modProfile.addFixedModification(getModification("TMT 10-plex of N-term"));
                prideParametersReport += "<br>" + "TMT 10-plex of N-term" + " (assumed fixed)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Phosphorylation")
                || pridePtmName.equalsIgnoreCase("L-aspartic 4-phosphoric anhydride")
                || pridePtmName.equalsIgnoreCase("O-phosphorylated residue")
                || pridePtmName.equalsIgnoreCase("Phospho")
                || pridePtmName.equalsIgnoreCase("phosphorylated residue")) {

            if (!modProfile.contains("Phosphorylation of S")) {
                modProfile.addVariableModification(getModification("Phosphorylation of S"));
                prideParametersReport += "<br>" + "Phosphorylation of S" + " (assumed variable)";
            }
            if (!modProfile.contains("Phosphorylation of T")) {
                modProfile.addVariableModification(getModification("Phosphorylation of T"));
                prideParametersReport += "<br>" + "Phosphorylation of T" + " (assumed variable)";
            }
            if (!modProfile.contains("Phosphorylation of Y")) {
                modProfile.addVariableModification(getModification("Phosphorylation of Y"));
                prideParametersReport += "<br>" + "Phosphorylation of Y" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Palmitoylation")) {

            if (!modProfile.contains("Palmitoylation of C")) {
                modProfile.addVariableModification(getModification("Palmitoylation of C"));
                prideParametersReport += "<br>" + "Palmitoylation of C" + " (assumed variable)";
            }
            if (!modProfile.contains("Palmitoylation of K")) {
                modProfile.addVariableModification(getModification("Palmitoylation of K"));
                prideParametersReport += "<br>" + "Palmitoylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Palmitoylation of S")) {
                modProfile.addVariableModification(getModification("Palmitoylation of S"));
                prideParametersReport += "<br>" + "Palmitoylation of S" + " (assumed variable)";
            }
            if (!modProfile.contains("Palmitoylation of T")) {
                modProfile.addVariableModification(getModification("Palmitoylation of T"));
                prideParametersReport += "<br>" + "Palmitoylation of T" + " (assumed variable)";
            }
            if (!modProfile.contains("Palmitoylation of protein N-term")) {
                modProfile.addVariableModification(getModification("Palmitoylation of protein N-term"));
                prideParametersReport += "<br>" + "Palmitoylation of protein N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Formylation")) {

            if (!modProfile.contains("Formylation of K")) {
                modProfile.addVariableModification(getModification("Formylation of K"));
                prideParametersReport += "<br>" + "Formylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Formylation of S")) {
                modProfile.addVariableModification(getModification("Formylation of S"));
                prideParametersReport += "<br>" + "Formylation of S" + " (assumed variable)";
            }
            if (!modProfile.contains("Formylation of T")) {
                modProfile.addVariableModification(getModification("Formylation of T"));
                prideParametersReport += "<br>" + "Formylation of T" + " (assumed variable)";
            }
            if (!modProfile.contains("Formylation of peptide N-term")) {
                modProfile.addVariableModification(getModification("Formylation of peptide N-term"));
                prideParametersReport += "<br>" + "Formylation of peptide N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Carbamylation")
                || pridePtmName.equalsIgnoreCase("carbamoylated residue")) {

            if (!modProfile.contains("Carbamylation of K")) {
                modProfile.addVariableModification(getModification("Carbamylation of K"));
                prideParametersReport += "<br>" + "Carbamylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Carbamilation of protein N-term")) {
                modProfile.addVariableModification(getModification("Carbamilation of protein N-term"));
                prideParametersReport += "<br>" + "Carbamilation of protein N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("3x(12)C labeled N6-propanoyl-L-lysine")) {

            if (!modProfile.contains("Propionyl of K light")) {
                modProfile.addVariableModification(getModification("Propionyl of K light"));
                prideParametersReport += "<br>" + "Propionyl of K light" + " (assumed variable)";
            }
            if (!modProfile.contains("Propionyl of peptide N-term light")) {
                modProfile.addVariableModification(getModification("Propionyl of peptide N-term light"));
                prideParametersReport += "<br>" + "Propionyl of peptide N-term light" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("3x(13)C labeled N6-propanoyl-L-lysine")) {

            if (!modProfile.contains("Propionyl of K heavy")) {
                modProfile.addVariableModification(getModification("Propionyl of K heavy"));
                prideParametersReport += "<br>" + "Propionyl of K heavy" + " (assumed variable)";
            }
            if (!modProfile.contains("Propionyl of peptide N-term heavy")) {
                modProfile.addVariableModification(getModification("Propionyl of peptide N-term heavy"));
                prideParametersReport += "<br>" + "Propionyl of peptide N-term heavy" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("3x(2)H residue methyl ester")) {

            if (!modProfile.contains("Trideuterated Methyl Ester of D")) {
                modProfile.addVariableModification(getModification("Trideuterated Methyl Ester of D"));
                prideParametersReport += "<br>" + "Trideuterated Methyl Ester of D" + " (assumed variable)";
            }
            if (!modProfile.contains("Trideuterated Methyl Ester of E")) {
                modProfile.addVariableModification(getModification("Trideuterated Methyl Ester of E"));
                prideParametersReport += "<br>" + "Trideuterated Methyl Ester of E" + " (assumed variable)";
            }
            if (!modProfile.contains("Trideuterated Methyl Ester of K")) {
                modProfile.addVariableModification(getModification("Trideuterated Methyl Ester of K"));
                prideParametersReport += "<br>" + "Trideuterated Methyl Ester of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Trideuterated Methyl Ester of R")) {
                modProfile.addVariableModification(getModification("Trideuterated Methyl Ester of R"));
                prideParametersReport += "<br>" + "Trideuterated Methyl Ester of R" + " (assumed variable)";
            }
            if (!modProfile.contains("Trideuterated Methyl Ester of peptide C-term")) {
                modProfile.addVariableModification(getModification("Trideuterated Methyl Ester of peptide C-term"));
                prideParametersReport += "<br>" + "Trideuterated Methyl Ester of peptide C-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("6x(13)C labeled residue")) {

            if (!modProfile.contains("Arginine 13C6")) {
                modProfile.addVariableModification(getModification("Arginine 13C6"));
                prideParametersReport += "<br>" + "Arginine 13C6" + " (assumed variable)";
            }
            if (!modProfile.contains("Lysine 13C6")) {
                modProfile.addVariableModification(getModification("Lysine 13C6"));
                prideParametersReport += "<br>" + "Lysine 13C6" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Acetyl")
                || pridePtmName.equalsIgnoreCase("N-acetylated residue")
                || pridePtmName.equalsIgnoreCase("N-acylated residue")
                || pridePtmName.equalsIgnoreCase("acetylated residue")) {

            if (!modProfile.contains("Acetylation of K")) {
                modProfile.addVariableModification(getModification("Acetylation of K"));
                prideParametersReport += "<br>" + "Acetylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Acetylation of peptide N-term")) {
                modProfile.addVariableModification(getModification("Acetylation of peptide N-term"));
                prideParametersReport += "<br>" + "Acetylation of peptide N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("dimethylated residue")) {

            if (!modProfile.contains("Dimethylation of K")) {
                modProfile.addVariableModification(getModification("Dimethylation of K"));
                prideParametersReport += "<br>" + "Dimethylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Dimethylation of R")) {
                modProfile.addVariableModification(getModification("Dimethylation of R"));
                prideParametersReport += "<br>" + "Dimethylation of R" + " (assumed variable)";
            }
            if (!modProfile.contains("Dimethylation of peptide N-term")) {
                modProfile.addVariableModification(getModification("Dimethylation of peptide N-term"));
                prideParametersReport += "<br>" + "Dimethylation of peptide N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("trimethylated residue")) {

            if (!modProfile.contains("Trimethylation of K")) {
                modProfile.addVariableModification(getModification("Trimethylation of K"));
                prideParametersReport += "<br>" + "Trimethylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Trimethylation of R")) {
                modProfile.addVariableModification(getModification("Trimethylation of R"));
                prideParametersReport += "<br>" + "Trimethylation of R" + " (assumed variable)";
            }
            if (!modProfile.contains("Trimethylation of protein N-term A")) {
                modProfile.addVariableModification(getModification("Trimethylation of protein N-term A"));
                prideParametersReport += "<br>" + "Trimethylation of protein N-term A" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Acetylation")) {

            if (!modProfile.contains("Acetylation of K")) {
                modProfile.addVariableModification(getModification("Acetylation of K"));
                prideParametersReport += "<br>" + "Acetylation of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Acetylation of peptide N-term")) {
                modProfile.addVariableModification(getModification("Acetylation of peptide N-term"));
                prideParametersReport += "<br>" + "Acetylation of peptide N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Guanidination")) {

            if (!modProfile.contains("Guanidination of K")) {
                modProfile.addVariableModification(getModification("Guanidination of K"));
                prideParametersReport += "<br>" + "Guanidination of K" + " (assumed variable)";
            }
            if (!modProfile.contains("Guanidination of peptide N-term")) {
                modProfile.addVariableModification(getModification("Guanidination of peptide N-term"));
                prideParametersReport += "<br>" + "Guanidination of peptide N-term" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Methylthio")) {

            if (!modProfile.contains("Methylthio of N")) {
                modProfile.addVariableModification(getModification("Methylthio of N"));
                prideParametersReport += "<br>" + "Methylthio of N" + " (assumed variable)";
            }
            if (!modProfile.contains("Methylthio of D")) {
                modProfile.addVariableModification(getModification("Methylthio of D"));
                prideParametersReport += "<br>" + "Methylthio of D" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Sulfo")
                || pridePtmName.equalsIgnoreCase("sulfated residue")) {

            if (!modProfile.contains("Sulfonation of Y")) {
                modProfile.addVariableModification(getModification("Sulfonation of Y"));
                prideParametersReport += "<br>" + "Sulfonation of Y" + " (assumed variable)";
            }
            if (!modProfile.contains("Sulfonation of S")) {
                modProfile.addVariableModification(getModification("Sulfonation of S"));
                prideParametersReport += "<br>" + "Sulfonation of S" + " (assumed variable)";
            }
            if (!modProfile.contains("Sulfonation of T")) {
                modProfile.addVariableModification(getModification("Sulfonation of T"));
                prideParametersReport += "<br>" + "Sulfonation of T" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Deamination")
                || pridePtmName.equalsIgnoreCase("Deamidated")
                || pridePtmName.equalsIgnoreCase("deamidated L-glutamine")
                || pridePtmName.equalsIgnoreCase("deamidated residue")
                || pridePtmName.equalsIgnoreCase("deaminated residue")) {

            if (!modProfile.contains("Deamidation of N")) {
                modProfile.addVariableModification(getModification("Deamidation of N"));
                prideParametersReport += "<br>" + "Deamidation of N" + " (assumed variable)";
            }
            if (!modProfile.contains("Deamidation of Q")) {
                modProfile.addVariableModification(getModification("Deamidation of Q"));
                prideParametersReport += "<br>" + "Deamidation of Q" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("Dioxidation")) {

            if (!modProfile.contains("Dioxidation of M")) {
                modProfile.addVariableModification(getModification("Dioxidation of M"));
                prideParametersReport += "<br>" + "Dioxidation of M" + " (assumed variable)";
            }
            if (!modProfile.contains("Dioxidation of W")) {
                modProfile.addVariableModification(getModification("Dioxidation of W"));
                prideParametersReport += "<br>" + "Dioxidation of W" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("dehydrated residue")
                || pridePtmName.equalsIgnoreCase("Dehydratation")) {

            if (!modProfile.contains("Dehydration of S")) {
                modProfile.addVariableModification(getModification("Dehydration of S"));
                prideParametersReport += "<br>" + "Dehydration of S" + " (assumed variable)";
            }
            if (!modProfile.contains("Dehydration of T")) {
                modProfile.addVariableModification(getModification("Dehydration of T"));
                prideParametersReport += "<br>" + "Dehydration of T" + " (assumed variable)";
            }

        } else if (pridePtmName.equalsIgnoreCase("No Modifications are included in the dataset")) {
            // ignore
        } else {

            // single modification mapping
            String utilitiesPtmName = getUtilitiesModificationName(pridePtmName);

            if (utilitiesPtmName != null) {
                if (!modProfile.contains(utilitiesPtmName)) {
                    if (isFixed) {
                        modProfile.addFixedModification(getModification(utilitiesPtmName));
                        prideParametersReport += "<br>" + utilitiesPtmName + " (assumed fixed)";
                    } else {
                        modProfile.addVariableModification(getModification(utilitiesPtmName));
                        prideParametersReport += "<br>" + utilitiesPtmName + " (assumed variable)";
                    }
                }
            } else if (!unknownPtms.contains(pridePtmName)) {
                unknownPtms.add(pridePtmName);
            }
        }

        return prideParametersReport;
    }

    /**
     * Tries to convert a modification name to a utilities modification
     * name. Null if no mapping is found.
     *
     * @param modName the modification name to map
     * @return the utilities modification name, or null if there is no mapping
     */
    public String getUtilitiesModificationName(String modName) {

        if (modName.equalsIgnoreCase("Carbamidomethyl")
                || modName.equalsIgnoreCase("S-carboxamidomethyl-L-cysteine")
                || modName.equalsIgnoreCase("iodoacetamide - site C")
                || modName.equalsIgnoreCase("iodoacetamide -site C")
                || modName.equalsIgnoreCase("iodoacetamide derivatized residue")
                || modName.equalsIgnoreCase("Iodoacetamide derivative")) {
            return "Carbamidomethylation of C";
        } else if (modName.equalsIgnoreCase("Oxidation")
                || modName.equalsIgnoreCase("monohydroxylated residue")
                || modName.equalsIgnoreCase("oxidized residue")) {
            return "Oxidation of M";
        } else if (modName.equalsIgnoreCase("Amidation")) {
            return "Amidation of peptide C-term";
        } else if (modName.equalsIgnoreCase("Carboxymethyl")
                || modName.equalsIgnoreCase("S-carboxymethyl-L-cysteine")
                || modName.equalsIgnoreCase("iodoacetic acid derivatized residue")) {
            return "Carboxymethylation of C";
        } else if (modName.equalsIgnoreCase("Farnesylation")) {
            return "Farnesylation of C";
        } else if (modName.equalsIgnoreCase("Geranyl-geranyl")) {
            return "Geranyl-geranyl of C";
        } else if (modName.equalsIgnoreCase("Homoserine")) {
            return "Homoserine of peptide C-term M";
        } else if (modName.equalsIgnoreCase("Homoserine lactone")) {
            return "Homoserine lactone of peptide C-term M";
        } else if (modName.equalsIgnoreCase("ICAT-C")
                || modName.equalsIgnoreCase("Applied Biosystems cleavable ICAT(TM) light")) {
            return "ICAT-O";
        } else if (modName.equalsIgnoreCase("ICAT-C:13C(9)")
                || modName.equalsIgnoreCase("Applied Biosystems cleavable ICAT(TM) heavy")) {
            return "ICAT-9";
        } else if (modName.equalsIgnoreCase("Lipoyl")) {
            return "Lipoyl of K";
        } else if (modName.equalsIgnoreCase("NIPCAM(C)")) {
            return "NIPCAM of C";
        } else if (modName.equalsIgnoreCase("Phosphopantetheine")) {
            return "phosphopantetheine s";
        } else if (modName.equalsIgnoreCase("Propionamide(C)")
                || modName.equalsIgnoreCase("Acrylamide adduct")) {
            return "Propionamide of C";
        } else if (modName.equalsIgnoreCase("Pyridylethyl")) {
            return "Pyridylethyl of C";
        } else if (modName.equalsIgnoreCase("(18)O label at both C-terminal oxygens")) {
            return "18O(2) of peptide C-term";
        } else if (modName.equalsIgnoreCase("(18)O monosubstituted residue")) {
            return "18O(1) of peptide C-term";
        } else if (modName.equalsIgnoreCase("(4,4,5,5-(2)H4)-L-lysine")) {
            return "Lysine 2H4";
        } else if (modName.equalsIgnoreCase("2-pyrrolidone-5-carboxylic acid (Gln)")
                || modName.equalsIgnoreCase("Ammonia-loss")) {
            return "Pyrolidone from Q";
        } else if (modName.equalsIgnoreCase("2-pyrrolidone-5-carboxylic acid (Glu)")
                || modName.equalsIgnoreCase("Glu->pyro-Glu")) {
            return "Pyrolidone from E";
        } else if (modName.equalsIgnoreCase("3-hydroxy-L-proline")) {
            return "Oxidation of P";
        } else if (modName.equalsIgnoreCase("3x(2)H labeled L-aspartic acid 4-methyl ester")) {
            return "Trideuterated Methyl Ester of D";
        } else if (modName.equalsIgnoreCase("4x(2)H labeled alpha-dimethylamino N-terminal residue")) {
            return "Dimethylation of peptide N-term 2H(4)";
        } else if (modName.equalsIgnoreCase("4x(2)H labeled dimethylated L-lysine")) {
            return "Dimethylation of K 2H(4)";
        } else if (modName.equalsIgnoreCase("5-methyl-L-arginine")) {
            return "Methylation of R";
        } else if (modName.equalsIgnoreCase("6x(13)C labeled L-arginine")) {
            return "Arginine 13C6";
        } else if (modName.equalsIgnoreCase("6x(13)C,4x(15)N labeled L-arginine")) {
            return "Arginine 13C6 15N4";
        } else if (modName.equalsIgnoreCase("6x(13)C labeled L-lysine")) {
            return "Lysine 13C6";
        } else if (modName.equalsIgnoreCase("6x(13)C,2x(15)N labeled L-lysine")) {
            return "Lysine 13C6 15N2";
        } else if (modName.equalsIgnoreCase("L-aspartic acid 4-methyl ester")) {
            return "Methylation of D";
        } else if (modName.equalsIgnoreCase("L-cysteic acid (L-cysteine sulfonic acid)")) {
            return "Oxidation of C";
        } else if (modName.equalsIgnoreCase("L-cysteine glutathione disulfide")) {
            return "Glutathione of C";
        } else if (modName.equalsIgnoreCase("L-cysteine methyl disulfide")
                || modName.equalsIgnoreCase("methyl methanethiosulfonate")) {
            return "Methylthio of C";
        } else if (modName.equalsIgnoreCase("L-cystine (cross-link)")) {
            return "Didehydro of T";
        } else if (modName.equalsIgnoreCase("L-glutamic acid 5-methyl ester (Glu)")
                || modName.equalsIgnoreCase("methylated glutamic acid")) {
            return "Methylation of E";
        } else if (modName.equalsIgnoreCase("L-homoarginine")) {
            return "Guanidination of K";
        } else if (modName.equalsIgnoreCase("L-methionine (R)-sulfoxide")
                || modName.equalsIgnoreCase("L-methionine (S)-sulfoxide")
                || modName.equalsIgnoreCase("L-methionine sulfoxide")) {
            return "Oxidation of M";
        } else if (modName.equalsIgnoreCase("L-methionine sulfone")) {
            return "Dioxidation of M";
        } else if (modName.equalsIgnoreCase("N-acetyl-L-asparagine")
                || modName.equalsIgnoreCase("N-acetyl-L-cysteine")
                || modName.equalsIgnoreCase("N-acetyl-L-glutamic acid")
                || modName.equalsIgnoreCase("N-acetyl-L-isoleucine")
                || modName.equalsIgnoreCase("N-acetyl-L-serine")
                || modName.equalsIgnoreCase("N-acetyl-L-tyrosine")
                || modName.equalsIgnoreCase("N2-acetyl-L-tryptophan")
                || modName.equalsIgnoreCase("alpha-amino acetylated residue")) {
            return "Acetylation of protein N-term";
        } else if (modName.equalsIgnoreCase("N-acetylated L-lysine")
                || modName.equalsIgnoreCase("N6-acetyl-L-lysine")) {
            return "Acetylation of K";
        } else if (modName.equalsIgnoreCase("N-ethylmaleimide derivatized cysteine")) {
            return "Nethylmaleimide of C";
        } else if (modName.equalsIgnoreCase("N-formyl-L-methionine")) {
            return "FormylMet of protein N-term";
        } else if (modName.equalsIgnoreCase("N-formylated residue")) {
            return "Formylation of peptide N-term"; // note: could also be the other formylations
        } else if (modName.equalsIgnoreCase("N-methyl-L-serine")) {
            return "Methylation of S";
        } else if (modName.equalsIgnoreCase("N6,N6-dimethyl-L-lysine")) {
            return "Dimethylation of K";
        } else if (modName.equalsIgnoreCase("N6-formyl-L-lysine")) {
            return "Formylation of K";
        } else if (modName.equalsIgnoreCase("N6-methyl-L-lysine")
                || modName.equalsIgnoreCase("methylated lysine")
                || modName.equalsIgnoreCase("monomethylated L-lysine")) {
            return "Methylation of K";
        } else if (modName.equalsIgnoreCase("N6-propanoyl-L-lysine")) {
            return "Propionyl of K light";
        } else if (modName.equalsIgnoreCase("O-(N-acetylamino)glucosyl-L-serine")) {
            return "HexNAc of S";
        } else if (modName.equalsIgnoreCase("O-(N-acetylamino)glucosyl-L-threonine")) {
            return "HexNAc of T";
        } else if (modName.equalsIgnoreCase("O-phospho-L-serine")) {
            return "Phosphorylation of S";
        } else if (modName.equalsIgnoreCase("O-phospho-L-threonine")) {
            return "Phosphorylation of T";
        } else if (modName.equalsIgnoreCase("O4&apos;-phospho-L-tyrosine")) {
            return "Phosphorylation of Y";
        } else if (modName.equalsIgnoreCase("S-carboxamidoethyl-L-cysteine")) {
            return "Propionamide of C";
        } else if (modName.equalsIgnoreCase("S-methyl-L-cysteine")) {
            return "Methylation of C";
        } else if (modName.equalsIgnoreCase("alpha-amino dimethylated residue")) {
            return "Dimethylation of N-term";
        } else if (modName.equalsIgnoreCase("amidated residue")) {
            return "Amidation of peptide C-term";
        } else if (modName.equalsIgnoreCase("deamidated L-asparagine")
                || modName.equalsIgnoreCase("deglycosylated asparagine")) {
            return "Deamidation of N";
        } else if (modName.equalsIgnoreCase("dihydroxylated residue - site W")) {
            return "Dioxidation of W";
        } else if (modName.equalsIgnoreCase("diiodinated residue")) {
            return "Diiodination of Y";
        } else if (modName.equalsIgnoreCase("hydroxylated lysine")
                || modName.equalsIgnoreCase("monohydroxylated lysine")) {
            return "Oxidation of K";
        } else if (modName.equalsIgnoreCase("iodoacetamide -site E")
                || modName.equalsIgnoreCase("iodoacetamide - site E")) {
            return "Carbamidomethylat\"iodoacetamide -site E\"ion of E";
        } else if (modName.equalsIgnoreCase("iodoacetamide N6-derivatized lysine")) {
            return "Carbamidomethylation of K";
        } else if (modName.equalsIgnoreCase("monomethylated L-aspartic acid")) {
            return "Methylation of D";
        } else if (modName.equalsIgnoreCase("thioacylation of primary amines - site N-term")) {
            return "Thioacyl of peptide N-term";
        } else if (modName.equalsIgnoreCase("ubiquitination signature dipeptidyl lysine")) {
            return "Ubiquitination of K";
        } else if (modName.equalsIgnoreCase("Label:13C(6)15N(2)")) {
            return "Lysine 13C(6) 15N(2)";
        } else {
            return null;
        }
    }

    /**
     * Returns the folder where the factory is saved.
     *
     * @return the folder where the factory is saved
     */
    public static String getSerializationFolder() {
        return SERIALIZATION_FILE_FOLDER;
    }

    /**
     * Sets the folder where the factory is saved.
     *
     * @param serializationFolder the folder where the factory is saved
     */
    public static void setSerializationFolder(String serializationFolder) {
        ModificationFactory.SERIALIZATION_FILE_FOLDER = serializationFolder;
    }

    /**
     * Sets the default modifications.
     */
    private void setDefaultModifications() {

        // Acetylation of K
        AtomChain atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        AtomChain atomChainRemoved = new AtomChain();
        AminoAcidPattern aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        String modificationName = "Acetylation of K";
        Modification modification = new Modification(
                ModificationType.modaa, modificationName, "ace",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1", "Acetyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00723", "N-acetylated L-lysine", null));
        modification.addReporterIon(ReporterIon.ACE_K_126);
        modification.addReporterIon(ReporterIon.ACE_K_143);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Acetylation of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Acetylation of peptide N-term"; // note: if name changed also change in TandemProcessBuilder
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "ace",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1", "Acetyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01458", "alpha-amino acetylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);        

        // Acetylation of protein N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Acetylation of protein N-term"; // note: if name changed also change in TandemProcessBuilder of SearchGUI and PsmImporter of PeptideShaker
        modification = new Modification(
                ModificationType.modn_protein, modificationName, "ace",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1", "Acetyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01458", "alpha-amino acetylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamidomethylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Carbamidomethylation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cmm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:4", "Carbamidomethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01060", "S-carboxamidomethyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamidomethylation of E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Carbamidomethylation of E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cmm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:4", "Carbamidomethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01216", "iodoacetamide derivatized glutamic acid", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamidomethylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Carbamidomethylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cmm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:4", "Carbamidomethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01212", "iodoacetamide N6-derivatized lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Oxidation of M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "Oxidation of M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:35", "Oxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00719", "L-methionine sulfoxide", null)); // @TODO: could also map to MOD:00425?
        modification.addNeutralLoss(NeutralLoss.CH4OS);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Oxidation of P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "Oxidation of P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:35", "Oxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00039", "4-hydroxy-L-proline", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Oxidation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Oxidation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:35", "Oxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01047", "monohydroxylated lysine", null)); // @TODO: maps to parent term "monohydroxylated lysine"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Oxidation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Oxidation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:35", "Oxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00210", "L-cysteine sulfenic acid", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Oxidation of N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Oxidation of N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:35", "Oxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01688", "3-hydroxy-L-asparagine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dioxydation of M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "Dioxidation of M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "diox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:425", "Dioxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00428", "dihydroxylated residue", null)); // @TODO: maps to parent term "dihydroxylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dioxydation of W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "Dioxidation of W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "diox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:425", "Dioxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00428", "dihydroxylated residue", null)); // @TODO: maps to parent term "dihydroxylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trioxidation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Trioxidation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "triox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:345", "Trioxidation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00460", "L-cysteic acid (L-cysteine sulfonic acid)", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Phosphorylation of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.P, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Phosphorylation of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "p",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:21", "Phospho", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00046", "O-phospho-L-serine", null));
        modification.addNeutralLoss(NeutralLoss.H3PO4);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Phosphorylation of T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.P, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Phosphorylation of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "p",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:21", "Phospho", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00047", "O-phospho-L-threonine", null));
        modification.addNeutralLoss(NeutralLoss.H3PO4);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Phosphorylation of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.P, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Phosphorylation of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "p",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:21", "Phospho", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00048", "O4'-phospho-L-tyrosine", null));
        modification.addNeutralLoss(NeutralLoss.HPO3);
        modification.addReporterIon(ReporterIon.PHOSPHO_Y);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Arg6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Arginine 13C(6)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01331", "6x(13)C labeled L-arginine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Arg10
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Arginine 13C(6) 15N(4)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:267", "Label:13C(6)15N(4)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00587", "6x(13)C,4x(15)N labeled L-arginine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Lys4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Lysine 2H(4)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:481", "Label:2H(4)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00942", "(4,4,5,5-(2)H4)-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Lys6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Lysine 13C(6)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01334", "6x(13)C labeled L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Lys8
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Lysine 13C(6) 15N(2)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:259", "Label:13C(6)15N(2)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00582", "6x(13)C,2x(15)N labeled L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Pro5
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 5);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "Proline 13C(5)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:772", "Label:13C(5)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01297", "5x(13)C labeled L-proline", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // 4-Hydroxyloproline
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 9);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "4-Hydroxyproline";
        modification = new Modification(
                ModificationType.modaa, modificationName, "hydroxy",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Leu7
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "Leucine 13C(6) 15N(1)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:695", "Label:13C(6)15N(1)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01285", "6x(13)C,1x(15)N labeled L-leucine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Ile7
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "Isoleucine 13C(6) 15N(1)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:695", "Label:13C(6)15N(1)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01286", "6x(13)C,1x(15)N labeled L-isoleucine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Label of K 2H4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Label of K 2H(4)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "2H(4)",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:481", "Label:2H(4)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00942", "(4,4,5,5-(2)H4)-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of K 2H4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Dimethylation of K 2H(4)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dimeth4",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:199 ", "Dimethyl:2H(4)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01254", "4x(2)H labeled dimethylated L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of K 2H6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Dimethylation of K 2H(6)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dimeth6",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1291", "Dimethyl:2H(6)", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of K 2H(6) 13C(2)
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Dimethylation of K 2H(6) 13C(2)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dimeth8",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:330", " Dimethyl:2H(6)13C(2)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of peptide N-term 2H4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Dimethylation of peptide N-term 2H(4)";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "dimeth4",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:199 ", "Dimethyl:2H(4)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of peptide N-term 2H6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Dimethylation of peptide N-term 2H(6)";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "dimeth6",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD: 1291", "Dimethyl:2H(6)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of peptide N-term 2H(6) 13C(2)
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Dimethylation of peptide N-term 2H(6) 13C(2)";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "dimeth8",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:330", " Dimethyl:2H(6)13C(2)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // 18O(2) of peptide C-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 2), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "18O(2) of peptide C-term";
        modification = new Modification(
                ModificationType.modc_peptide, modificationName, "18O",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD: 193", "Label:18O(2)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00546", "(18)O label at both C-terminal oxygens", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // 18O(1) of peptide C-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 2), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 1);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "18O(1) of peptide C-term";
        modification = new Modification(
                ModificationType.modc_peptide, modificationName, "18O",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:258", "Label:18O(1)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00581", "(18)O monosubstituted residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICAT-0
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 17);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "ICAT-O";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:105", "ICAT-C", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00480", "Applied Biosystems cleavable ICAT(TM) light", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICAT-9
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 9);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 17);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "ICAT-9";
        modification = new Modification(
                ModificationType.modaa, modificationName, "*",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:106", "ICAT-C:13C(9)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00481", "Applied Biosystems cleavable ICAT(TM) heavy", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL0 of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "ICPL0 of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "icpl0",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:365", "ICPL", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01230", "Bruker Daltonics SERVA-ICPL(TM) quantification chemistry, light form - site K", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL0 of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "ICPL0 of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "icpl0",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:365", "ICPL", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL4 of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "ICPL4 of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "icpl4",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:687", "ICPL:2H(4)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01359", "Bruker Daltonics SERVA-ICPL(TM) quantification chemistry, medium form - site K", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL4 of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "ICPL4 of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "icpl4",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:687", "ICPL:2H(4)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01358", "Bruker Daltonics SERVA-ICPL(TM) quantification chemistry, medium form - site N-term", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL6 of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "ICPL6 of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "icpl6",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:364", "ICPL:13C(6)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL6 of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "ICPL6 of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "icpl6",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:364", "ICPL:13C(6)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL10 of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "ICPL10 of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "icpl10",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:866", "ICPL:13C(6)2H(4)", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01287", "Bruker Daltonics SERVA-ICPL(TM) quantification chemistry, heavy form - site K", null)); // @TODO: the mass in Unimod and PSI-MOD is not the same!
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ICPL10 of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "ICPL10 of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "icpl10",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:866", "ICPL:13C(6)2H(4)", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // mTRAQ of K light
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "mTRAQ of K light";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mTRAQ0",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD: 888", "mTRAQ", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01863", "mTRAQ reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "mTRAQ reporter+balance reagent acylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // mTRAQ of peptide N-term light
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "mTRAQ of peptide N-term light";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "mTRAQ0",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:888", "mTRAQ", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01863", "mTRAQ reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "mTRAQ reporter+balance reagent acylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // mTRAQ of K 13C3 15N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "mTRAQ of K 13C(3) 15N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mTRAQ4",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:889", "mTRAQ:13C(3)15N(1) ", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01863", "mTRAQ reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "mTRAQ reporter+balance reagent acylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // mTRAQ of peptide N-term 13C3 15N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "mTRAQ of peptide N-term 13C(3) 15N";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "mTRAQ4",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:889", "mTRAQ:13C(3)15N(1) ", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01863", "mTRAQ reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "mTRAQ reporter+balance reagent acylated residue"...

        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // mTRAQ of K 13C6 15N2
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "mTRAQ of 13C(6) 15N(2)";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mTRAQ8",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1302", "mTRAQ:13C(6)15N(2) ", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01863", "mTRAQ reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "mTRAQ reporter+balance reagent acylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // mTRAQ of peptide N-term 13C3 15N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "mTRAQ of peptide N-term 13C(6) 15N(2)";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "mTRAQ8",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1302", "mTRAQ:13C(6)15N(2) ", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01863", "mTRAQ reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "mTRAQ reporter+balance reagent acylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iTRAQ 4-plex of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "iTRAQ 4-plex of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "iTRAQ",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:214", "iTRAQ4plex", null)); // @TODO: check cv term and mass!!!
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01518", "iTRAQ4plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "iTRAQ4plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_114);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_115);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_116);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_117);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iTRAQ 4-plex of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "iTRAQ 4-plex of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "iTRAQ",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:214", "iTRAQ4plex", null)); // @TODO: check cv term and mass!!!
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01518", "iTRAQ4plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "iTRAQ4plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_114);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_115);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_116);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_117);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iTRAQ 4-plex of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "iTRAQ 4-plex of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "iTRAQ",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:214", "iTRAQ4plex", null)); // @TODO: check cv term and mass!!!
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01518", "iTRAQ4plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "iTRAQ4plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_114);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_115);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_116);
        modification.addReporterIon(ReporterIon.iTRAQ4Plex_117);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iTRAQ 8-plex of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 24);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "iTRAQ 8-plex of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "iTRAQ",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:730", "iTRAQ8plex", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01526", "iTRAQ8plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "iTRAQ8plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_113);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_114);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_115);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_116);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_117);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_118);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_119);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_121);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iTRAQ 8-plex of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 24);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "iTRAQ 8-plex of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "iTRAQ",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:730", "iTRAQ8plex", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01526", "iTRAQ8plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "iTRAQ8plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_113);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_114);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_115);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_116);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_117);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_118);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_119);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_121);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iTRAQ 8-plex of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 24);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "iTRAQ 8-plex of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "iTRAQ",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:730", "iTRAQ8plex", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01526", "iTRAQ8plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "iTRAQ8plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_113);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_114);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_115);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_116);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_117);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_118);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_119);
        modification.addReporterIon(ReporterIon.iTRAQ8Plex_121);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 2-plex of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 11);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "TMT 2-plex of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:738", "TMT2plex", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 2-plex of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 11);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 2-plex of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:738", "TMT2plex", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        /// TMT 6-plex of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = null;
        modificationName = "TMT 6-plex of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:737", "TMT6plex", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01715", "TMT6plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "TMT6plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 6-plex of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 6-plex of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:737", "TMT6plex", null)); // note: does not have a PSI name, using interim name
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01715", "TMT6plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "TMT6plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 6-plex + K+4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 6-plex of K+4";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 6-plex of K+6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 6-plex of K+6";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 6-plex of K+8
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 6-plex of K+8";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 10-plex of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = null;
        modificationName = "TMT 10-plex of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:737", "TMT6plex", null)); // note: "PSI: Alt. Description: Also applies to TMT10plex", i.e., no term for 10-plex
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01715", "TMT6plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "TMT6plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 10-plex of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 10-plex of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:737", "TMT6plex", null)); // note: "PSI: Alt. Description: Also applies to TMT10plex", i.e., no term for 10-plex
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 10-plex + K+4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 10-plex of K+4";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        /// @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 10-plex of K+6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 10-plex of K+6";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 10-plex of K+8
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 10-plex of K+8";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        //// @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 11-plex of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = null;
        modificationName = "TMT 11-plex of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:737", "TMT6plex", null)); // note: "PSI: Alt. Description: Also applies to TMT10plex", i.e., no term for 10-plex (no mention of 11-plex though...)
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01715", "TMT6plex reporter+balance reagent acylated residue", null)); // @TODO: maps to parent term "TMT6plex reporter+balance reagent acylated residue"...
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_131C);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        modification.addReporterIon(ReporterIon.TMT_131C_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 11-plex of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 11-plex of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:737", "TMT6plex", null)); // note: "PSI: Alt. Description: Also applies to TMT10plex", i.e., no term for 10-plex (no mention of 11-plex though...)
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_131C);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        modification.addReporterIon(ReporterIon.TMT_131C_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 11-plex + K+4
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 11-plex of K+4";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        /// @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        modification.addReporterIon(ReporterIon.TMT_131C_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 11-plex of K+6
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 11-plex of K+6";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        modification.addReporterIon(ReporterIon.TMT_131C_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMT 11-plex of K+8
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 20);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMT 11-plex of K+8";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        modification.addReporterIon(ReporterIon.TMT_126);
        modification.addReporterIon(ReporterIon.TMT_127C);
        modification.addReporterIon(ReporterIon.TMT_127N);
        modification.addReporterIon(ReporterIon.TMT_128C);
        modification.addReporterIon(ReporterIon.TMT_128N);
        modification.addReporterIon(ReporterIon.TMT_129C);
        modification.addReporterIon(ReporterIon.TMT_129N);
        modification.addReporterIon(ReporterIon.TMT_130C);
        modification.addReporterIon(ReporterIon.TMT_130N);
        modification.addReporterIon(ReporterIon.TMT_131);
        modification.addReporterIon(ReporterIon.TMT_126_ETD);
        modification.addReporterIon(ReporterIon.TMT_127C_ETD);
        modification.addReporterIon(ReporterIon.TMT_127N_ETD);
        modification.addReporterIon(ReporterIon.TMT_128C_ETD);
        modification.addReporterIon(ReporterIon.TMT_128N_ETD);
        modification.addReporterIon(ReporterIon.TMT_129C_ETD);
        modification.addReporterIon(ReporterIon.TMT_129N_ETD);
        modification.addReporterIon(ReporterIon.TMT_130C_ETD);
        modification.addReporterIon(ReporterIon.TMT_130N_ETD);
        modification.addReporterIon(ReporterIon.TMT_131_ETD);
        modification.addReporterIon(ReporterIon.TMT_131C_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMTpro of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 25);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = null;
        modificationName = "TMTpro of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:2016", "TMTpro", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI_MOD cv term?
        modification.addReporterIon(ReporterIon.TMTpro_126);
        modification.addReporterIon(ReporterIon.TMTpro_127N);
        modification.addReporterIon(ReporterIon.TMTpro_127C);
        modification.addReporterIon(ReporterIon.TMTpro_128N);
        modification.addReporterIon(ReporterIon.TMTpro_128C);
        modification.addReporterIon(ReporterIon.TMTpro_129N);
        modification.addReporterIon(ReporterIon.TMTpro_129C);
        modification.addReporterIon(ReporterIon.TMTpro_130N);
        modification.addReporterIon(ReporterIon.TMTpro_130C);
        modification.addReporterIon(ReporterIon.TMTpro_131N);
        modification.addReporterIon(ReporterIon.TMTpro_131C);
        modification.addReporterIon(ReporterIon.TMTpro_132N);
        modification.addReporterIon(ReporterIon.TMTpro_132C);
        modification.addReporterIon(ReporterIon.TMTpro_133N);
        modification.addReporterIon(ReporterIon.TMTpro_133C);
        modification.addReporterIon(ReporterIon.TMTpro_134N);
        modification.addReporterIon(ReporterIon.TMTpro_134C);
        modification.addReporterIon(ReporterIon.TMTpro_135N);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // TMTpro of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 25);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "TMTpro of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "TMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:2016", "TMTpro", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI_MOD cv term?
        modification.addReporterIon(ReporterIon.TMTpro_126);
        modification.addReporterIon(ReporterIon.TMTpro_127N);
        modification.addReporterIon(ReporterIon.TMTpro_127C);
        modification.addReporterIon(ReporterIon.TMTpro_128N);
        modification.addReporterIon(ReporterIon.TMTpro_128C);
        modification.addReporterIon(ReporterIon.TMTpro_129N);
        modification.addReporterIon(ReporterIon.TMTpro_129C);
        modification.addReporterIon(ReporterIon.TMTpro_130N);
        modification.addReporterIon(ReporterIon.TMTpro_130C);
        modification.addReporterIon(ReporterIon.TMTpro_131N);
        modification.addReporterIon(ReporterIon.TMTpro_131C);
        modification.addReporterIon(ReporterIon.TMTpro_132N);
        modification.addReporterIon(ReporterIon.TMTpro_132C);
        modification.addReporterIon(ReporterIon.TMTpro_133N);
        modification.addReporterIon(ReporterIon.TMTpro_133C);
        modification.addReporterIon(ReporterIon.TMTpro_134N);
        modification.addReporterIon(ReporterIon.TMTpro_134C);
        modification.addReporterIon(ReporterIon.TMTpro_135N);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iodoTMT zero of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 28);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "iodoTMT zero of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "iodoTMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1341", "iodoTMT", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI_MOD cv term?
        modification.addReporterIon(ReporterIon.iodoTMT_zero);
        modification.addReporterIon(ReporterIon.iodoTMT_zero_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // iodoTMT 6-plex of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 28);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 1), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "iodoTMT 6-plex of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "iodoTMT",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1342", "iodoTMT", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI_MOD cv term?
        modification.addReporterIon(ReporterIon.iodoTMT_126);
        modification.addReporterIon(ReporterIon.iodoTMT_127);
        modification.addReporterIon(ReporterIon.iodoTMT_128);
        modification.addReporterIon(ReporterIon.iodoTMT_129);
        modification.addReporterIon(ReporterIon.iodoTMT_130);
        modification.addReporterIon(ReporterIon.iodoTMT_131);
        modification.addReporterIon(ReporterIon.iodoTMT_126_ETD);
        modification.addReporterIon(ReporterIon.iodoTMT_127_ETD);
        modification.addReporterIon(ReporterIon.iodoTMT_128_ETD);
        modification.addReporterIon(ReporterIon.iodoTMT_129_ETD);
        modification.addReporterIon(ReporterIon.iodoTMT_130_ETD);
        modification.addReporterIon(ReporterIon.iodoTMT_131_ETD);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Ubiquitination of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Ubiquitination of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:121", "GG", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00492", "ubiquitination signature dipeptidyl lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Methylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "meth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:34", "Methyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00085", "N6-methyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylation of R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Methylation of R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "meth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:34", "Methyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00414", "monomethylated L-arginine", null)); // @TODO: maps to parent term "monomethylated L-arginine"...
        modification.addReporterIon(ReporterIon.METHYL_R_70);
        modification.addReporterIon(ReporterIon.METHYL_R_87);
        modification.addReporterIon(ReporterIon.METHYL_R_112);
        modification.addReporterIon(ReporterIon.METHYL_R_115);
        modification.addReporterIon(ReporterIon.METHYL_R_143);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylation of E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Methylation of E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "meth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:34", "Methyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00081", "L-glutamic acid 5-methyl ester (Glu)", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Methylation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "meth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:34", "Methyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00660", "methylated cysteine", null)); // @TODO: maps to parent term "methylated cysteine"
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylation of D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Methylation of D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "meth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:34", "Methyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00079", "N4-methyl-L-asparagine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylation of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Methylation of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "meth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:34", "Methyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01782", "N-methyl-L-serine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Dimethylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:36", "Dimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00084", "N6,N6-dimethyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Dimethylation of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "dimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:36", "Dimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01686", "alpha-amino dimethylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Dimethylation of R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:36", "Dimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00783", "dimethylated L-arginine", null));
        modification.addReporterIon(ReporterIon.DI_METHYL_R_112);
        modification.addReporterIon(ReporterIon.DI_METHYL_R_115);
        modification.addReporterIon(ReporterIon.DI_METHYL_R_157);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dimethylation of N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Dimethylation of N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:36", "Dimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00316", "N4,N4-dimethyl-L-asparagine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trimethylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Trimethylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "trimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:37", "Trimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00855", "N6,N6,N6-trimethyl-L-lysine (from L-lysinium residue)", null));
        modification.addNeutralLoss(NeutralLoss.C3H9N);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trimethylation of R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Trimethylation of R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "trimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:37", "Trimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01669", "trimethyl-L-arginine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trimethylation of protein N-term A
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "Trimethylation of protein N-term A";
        modification = new Modification(
                ModificationType.modn_protein, modificationName, "trimeth",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:37", "Trimethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01687", "alpha-amino trimethylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Pyrolidone from E
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Pyrolidone from E"; // note: if name changed also change in TandemProcessBuilder of SearchGUI and PsmImporter of PeptideShaker
        modification = new Modification(
                ModificationType.modnaa_peptide, modificationName, "pyro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:27", "Glu->pyro-Glu", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00420", "2-pyrrolidone-5-carboxylic acid (Glu)", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Pyrolidone from Q
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Pyrolidone from Q"; // note: if name changed also change in TandemProcessBuilder of SearchGUI and PsmImporter of PeptideShaker
        modification = new Modification(
                ModificationType.modnaa_peptide, modificationName, "pyro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:28", "Gln->pyro-Glu", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00040", "2-pyrrolidone-5-carboxylic acid (Gln)", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Pyrolidone from carbamidomethylated C
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Pyrolidone from carbamidomethylated C"; // note: if name changed also change in TandemProcessBuilder of SearchGUI and PsmImporter of PeptideShaker
        modification = new Modification(
                ModificationType.modnaa_peptide, modificationName, "pyro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:385", "Ammonia-loss", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01871", "cyclized N-terminal S-carboxamidomethyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // HexNAc of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 13);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 5);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "HexNAc of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glyco",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:43", "HexNAc", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // HexNAc of T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 13);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 5);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "HexNAc of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glyco",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:43", "HexNAc", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Hex(1)NAc(1) of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 14);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 23);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 10);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Hex(1)NAc(1) of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glyco",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:793", "Hex(1)HexNAc(1) ", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Hex(1)NAc(1) of T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 14);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 23);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 10);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Hex(1)NAc(1) of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glyco",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:793", "Hex(1)HexNAc(1) ", null)); // note: does not have a PSI name, using interim name
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Hexose of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 10);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 5);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Hexose of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "hex",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:41", "Hex", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01347", "hexose glycated L-lysine", null));

        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Hex(5) HexNAc(4) NeuAc(2) of N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 84);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 136);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 61);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Hex(5) HexNAc(4) NeuAc(2) of N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glyco",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        // @TODO: add Unimod mapping?
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Hex(5) HexNAc(4) NeuAc(2) Na of N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 84);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 135);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 61);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.Na, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Hex(5) HexNAc(4) NeuAc(2) Na of N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glyco",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Glyco);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1408", "Hex(5)HexNAc(4)NeuAc(2) ", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // SUMO-2/3 Q87R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 18);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 29);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 8);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "SUMO-2/3 Q87R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sumo",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        // @TODO: add Unimod mapping?
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01149", "sumoylated lysine", null)); // @TODO: maps to parent term "sumoylated lysine"...
        modification.addReporterIon(ReporterIon.QQ);
        modification.addReporterIon(ReporterIon.QQ_H2O);
        modification.addReporterIon(ReporterIon.QQT);
        modification.addReporterIon(ReporterIon.QQT_H2O);
        modification.addReporterIon(ReporterIon.QQTG);
        modification.addReporterIon(ReporterIon.QQTG_H2O);
        modification.addReporterIon(ReporterIon.QQTGG);
        modification.addReporterIon(ReporterIon.QQTGG_H2O);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Deamidation of N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Deamidation of N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "deam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:7", "Deamidated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00684", "deamidated L-asparagine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Deamidation of Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Deamidation of Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "deam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:7", "Deamidated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00685", "deamidated L-glutamine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Deamidation of N 18O
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 2), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Deamidation of N 18O";
        modification = new Modification(
                ModificationType.modaa, modificationName, "deam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD: 366", "Deamidation in presence of O18", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00791", "1x(18)O labeled deamidated L-glutamine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamylation of protein N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Carbamilation of protein N-term";
        modification = new Modification(
                ModificationType.modn_protein, modificationName, "cm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:5", "Carbamyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01679", "alpha-aminocarbamoylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Carbamilation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:5", "Carbamyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01678", "N6-carbamoyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamylation of R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Carbamilation of R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:5", "Carbamyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00398", "carbamoylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Carbamilation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:5", "Carbamyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00337", "S-carbamoyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carbamylation of M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "Carbamilation of M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:5", "Carbamyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00398", "carbamoylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Acetaldehyde +26
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Acetaldehyde +26";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "ace",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Other);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:254", "Delta:H(2)C(2)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00577", "acetaldehyde +26", null)); // @TODO: this PSI-MOD not peptide n-term specific
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Sodium adduct to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Na, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Sodium adduct to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "Na",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:30", "Cation:Na", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01295", "monosodium L-aspartate", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Sodium adduct to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Na, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Sodium adduct to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "Na",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:30", "Cation:Na", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01296", "monosodium L-glutamate", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Amidation of peptide C-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 1);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Amidation of peptide C-term";
        modification = new Modification(
                ModificationType.modc_peptide, modificationName, "am",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:2", "Amidated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00883", "C1-amidated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Amidation of protein C-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 1);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Amidation of protein C-term";
        modification = new Modification(
                ModificationType.modc_protein, modificationName, "am",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:2", "Amidated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00883", "C1-amidated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Sulfation of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Sulfation of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "s",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:40", "Sulfo", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00366", "O-sulfo-L-serine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Sulfation of T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Sulfation of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "s",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:40", "Sulfo", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00367", "O-sulfo-L-threonine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Sulfation of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Sulfation of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "s",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:40", "Sulfo", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00181", "O4'-sulfo-L-tyrosine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Palmitoylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 30);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Palmitoylation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "palm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:47", "Palmitoylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00440", "palmitoylated residue", null)); // @TODO: maps to parent term "palmitoylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Palmitoylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 30);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Palmitoylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "palm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:47", "Palmitoylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00086", "N6-palmitoyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Palmitoylation of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 30);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Palmitoylation of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "palm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:47", "Palmitoylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00089", "O-palmitoyl-L-serine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Palmitoylation of T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 30);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Palmitoylation of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "palm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:47", "Palmitoylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00088", "O-palmitoyl-L-threonine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Palmitoylation of protein N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 30);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 16);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Palmitoylation of protein N-term";
        modification = new Modification(
                ModificationType.modn_protein, modificationName, "palm",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:47", "Palmitoylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01685", "alpha-amino palmitoylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Formylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Formylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "form",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:122", "Formylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00216", "N6-formyl-L-lysine", null));
        modification.addReporterIon(ReporterIon.FORMYL_K);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Formylation of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Formylation of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "form",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:122", "Formylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01222", "O-formyl-L-serine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Formylation of T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Formylation of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "form",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:122", "Formylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01221", "O-formyl-L-threonine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Formylation of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Formylation of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "form",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:122", "Formylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00409", "N-formylated residue", null)); // @TODO: maps to parent term "N-formylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Formylation of protein N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Formylation of protein N-term";
        modification = new Modification(
                ModificationType.modn_protein, modificationName, "form",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:122", "Formylation", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00409", "N-formylated residue", null)); // @TODO: maps to parent term "N-formylated residue"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionyl of K light
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Propionyl of K light";
        modification = new Modification(
                ModificationType.modaa, modificationName, "prop",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:58", "Propionyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01232", "3x(12)C labeled N6-propanoyl-L-lysine", null));

        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionyl of peptide N-term light
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Propionyl of peptide N-term light";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "prop",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:58", "Propionyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01895", "alpha-amino 3x(12)C-labeled propanoylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionyl of K heavy
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Propionyl of K heavy";
        modification = new Modification(
                ModificationType.modaa, modificationName, "prop",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:59", "Propionyl:13C(3)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01231", "3x(13)C labeled N6-propanoyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionyl of peptide N-term heavy
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Propionyl of peptide N-term heavy";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "prop",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:59", "Propionyl:13C(3)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00452", "alpha-amino 3x(13)C-labeled propanoylated residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trideuterated Methyl Ester of D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Trideuterated Methyl Ester of D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "methyl(d3)",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:298", "Methyl:2H(3)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01241", "3x(2)H labeled L-aspartic acid 4-methyl ester", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trideuterated Methyl Ester of E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Trideuterated Methyl Ester of E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "methyl(d3)",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:298", "Methyl:2H(3)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01242", "3x(2)H labeled L-glutamic acid 5-methyl ester", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trideuterated Methyl Ester of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Trideuterated Methyl Ester of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "methyl(d3)",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:298", "Methyl:2H(3)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trideuterated Methyl Ester of R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Trideuterated Methyl Ester of R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "methyl(d3)",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:298", "Methyl:2H(3)", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Trideuterated Methyl Ester of peptide C-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 1), 3);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Trideuterated Methyl Ester of peptide C-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "methyl(d3)",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Labeling);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:298", "Methyl:2H(3)", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00617", "3x(2)H residue methyl ester", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxymethylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Carboxymethylation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "carbox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:6", "Carboxymethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01061", "S-carboxymethyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxymethylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Carboxymethylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "carbox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:6", "Carboxymethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01094", "N6-carboxymethyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxymethylation of W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "Carboxymethylation of W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "carbox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:6", "Carboxymethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01328", "iodoacetic acid - site W", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxymethylation of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Carboxymethylation of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "carbox",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:6", "Carboxymethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01084", "iodoacetic acid derivatized amino-terminal residue", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Farnesylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 24);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 15);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Farnesylation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "far",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:44", "Farnesyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00111", "S-farnesyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Geranyl-geranyl of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 32);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 20);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Geranyl-geranyl of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "geranyl",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:48", "GeranylGeranyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00113", "S-geranylgeranyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Guanidination of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Guanidination of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "guan",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:52", "Guanidinyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00445", "L-homoarginine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Guanidination of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Guanidination of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "guan",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:52", "Guanidinyl", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Homoserine of peptide C-term M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "Homoserine of peptide C-term M";
        modification = new Modification(
                ModificationType.modcaa_peptide, modificationName, "hse",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:10", "Met->Hse", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00403", "homoserine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Homoserine lactone of peptide C-term M
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "Homoserine lactone of peptide C-term M";
        modification = new Modification(
                ModificationType.modcaa_peptide, modificationName, "hsel",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:11", "Met->Hsl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00404", "homoserine lactone", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Lipoyl of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 12);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Lipoyl of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "lip",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:42", "Lipoyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00127", "N6-lipoyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylthio of D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Methylthio of D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mmts",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:39", "Methylthio", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00237", "L-beta-methylthioaspartic acid", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylthio of N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Methylthio of N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mmts",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:39", "Methylthio", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00325", "L-beta-methylthioasparagine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Methylthio of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Methylthio of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mmts",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:39", "Methylthio", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00110", "L-cysteine methyl disulfide", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // NIPCAM of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 9);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "NIPCAM of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "nipcam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:17", "NIPCAM", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00410", "S-(N-isopropylcarboxamidomethyl)-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionamide of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Propionamide of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "propam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:24", "Propionamide", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00417", "S-carboxamidoethyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionamide of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Propionamide of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "propam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:24", "Propionamide", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Propionamide of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Propionamide of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "propam",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:24", "Propionamide", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Pyridylethyl of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Pyridylethyl of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "pyri",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:31", "Pyridylethyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00424", "S-pyridylethyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dehydration of S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "Dehydration of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dehyd",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:23", "Dehydrated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00189", "dehydroalanine (Ser)", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Dehydration of T
        atomChainAdded = new AtomChain();;
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Dehydration of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "dehyd",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:23", "Dehydrated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00190", "dehydrobutyrine (Thr)", null)); // @TODO: maps to parent term "dehydrobutyrine (Thr)"...
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Nethylmaleimide of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Nethylmaleimide of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "nem",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:108", "Nethylmaleimide", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00483", "N-ethylmaleimide derivatized cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Glutathione of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 15);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 10);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Glutathione of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glut",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:55", "Glutathione", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00234", "L-cysteine glutathione disulfide", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // FormylMet of protein N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 9);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "FormylMet of protein N-term";
        modification = new Modification(
                ModificationType.modn_protein, modificationName, "nmet",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:107", "FormylMet", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Didehydro of T
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "Didehydro of T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "didehyro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:401", "Didehydro", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01345", "2-amino-3-oxobutanoic acid", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Didehydro of Y
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Didehydro of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "didehyro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:401", "Didehydro", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00706", "dehydrogenated tyrosine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Thioacyl of peptide N-term
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.S, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = new AminoAcidPattern();
        modificationName = "Thioacyl of peptide N-term";
        modification = new Modification(
                ModificationType.modn_peptide, modificationName, "thioacyl",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:126", "Thioacyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01223", "thioacylation of primary amines - site N-term", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Diiodination of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.I, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Diiodination of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "diiodo",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:130", "Diiodo", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01140", "diiodinated tyrosine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Citrullination of R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "Citrullination of R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cit",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:7", "Deamidated", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00400", "deamidated residue", null));
        modification.addNeutralLoss(NeutralLoss.HCNO);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S-nitrosylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 1);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 1);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 1);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "S-nitrosylation";
        modification = new Modification(
                ModificationType.modaa, modificationName, "nitrosyl",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:275", "Nitrosyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00235", "S-nitrosyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Heme B of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 32);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 34);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.Fe, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Heme B of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "heme",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:390", "Heme", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Heme B of H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 32);
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 34);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.Fe, 0), 1);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "Heme B of H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "heme",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:390", "Heme", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Carboxylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "carb",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:299", "Carboxy", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00123", "N6-carboxy-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxylation of D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Carboxylation of D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "carb",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:299", "Carboxy", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00309", "L-beta-carboxyaspartic acid", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Carboxylation of E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Carboxylation of E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "carb",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:299", "Carboxy", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00041", "L-gamma-carboxyglutamic acid", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Nitrosylation of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Nitrosylation of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "nitro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:354", "Nitro", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01352", "nitrated L-tyrosine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Nitrosylation of C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "Nitrosylation of C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "nitro",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:275", "Nitrosyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00235", "S-nitrosyl-L-cysteine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Butyrylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Butyrylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "buty",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1289", "Butyryl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01781", "N6-butanoyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Crotonylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Crotonylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "croto",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1363", "Crotonyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01892", "N6-crotonyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Glutarylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Glutarylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "glur",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1848", "Gluratylation", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Malonylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Malonylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "malo",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:747", "Malonyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01893", "N6-malonyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Succinylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Succinylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "suc",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:64", "Succinyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01819", "N6-succinyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // ADP-ribosylation of S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 15);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 21);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 13);
        atomChainAdded.append(new AtomImpl(Atom.P, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "ADP-ribosylation of S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "adp",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:213", "ADP-Ribosyl", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00242", "O-(ADP-ribosyl)-L-serine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Pyridoxal phosphate of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.P, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Pyridoxal phosphate of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "pyri",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Biological);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:46", "PyridoxalPhosphate", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00128", "N6-pyridoxal phosphate-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Biotinylation of K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 10);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 14);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "Biotinylation of K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "biot",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:3", "Biotin", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00126", "N6-biotinyl-L-lysine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Oxidation to Kynurenine of W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "Oxidation to Kynurenine of W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "kynu",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:351", "kynurenin", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:00462", "L-kynurenine", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Quinone of Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Quinone of Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "quin",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:392", "Quinone", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Quinone of W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "Quinone of W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "quin",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Less_Common);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:392", "Quinone", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Potassium on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.K, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Potassium on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "pot",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:530", "Cation:K", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01245", "potassium L-aspartate", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Potassium on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.K, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Potassium on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "pot",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:530", "Cation:K", null));
        modification.setPsiModCvTerm(new CvTerm("MOD", "MOD:01244", "potassium L-glutamate", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Calcium on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Ca, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Calcium on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ca",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:951", "Cation:Ca[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Calcium on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Ca, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Calcium on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "ca",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:951", "Cation:Ca[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Zinc on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Zn, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Zinc on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "zn",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:954", "Cation:Zn[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Zinc on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Zn, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Zinc on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "zn",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:954", "Cation:Zn[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Iron[II] on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Fe, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Iron[II] on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "fe",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:952", "Cation:Fe[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Iron[II] on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Fe, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Iron[II] on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "fe",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:952", "Cation:Fe[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Iron[III] on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Fe, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Iron[III] on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "fe",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1870", "Cation:Fe[III]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Iron[III] on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Fe, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Iron[III] on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "fe",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1870", "Cation:Fe[III]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Magnesium on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Mg, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Magnesium on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mg",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:956", "Cation:Mg[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Magnesium on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Mg, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Magnesium on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "mg",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:956", "Cation:Mg[II]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Copper on D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Cu, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "Copper on D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cu",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:531", "Cation:Cu[I]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Copper on E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.Cu, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "Copper on E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "cu",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Metal);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:531", "Cation:Cu[I]", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Ammonia loss from N
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 1);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "Ammonia loss from N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "-nh3",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Common_Artifact);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:385", "Ammonia-loss", null));
        // @TODO: add PSI-MOD mapping?
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        addSubstitutions();
    }

    /**
     * Add the substitutions. Has to be done in a separate method or the
     * combined method becomes too large...
     */
    private void addSubstitutions() {

        // W to G
        AtomChain atomChainAdded = new AtomChain();
        AtomChain atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 9);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        AminoAcidPattern aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        String modificationName = "W to G";
        Modification modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:676", "Trp->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 8);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1224", "Trp->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1239", "Tyr->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 9);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:646", "Arg->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 8);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:673", "Trp->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1237", "Tyr->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1093", "Phe->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1232", "Trp->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:677", "Trp->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1189", "Arg->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1234", "Trp->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 8);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:674", "Trp->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1117", "His->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1090", "Phe->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:679", "Tyr->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1147", "Met->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1231", "Trp->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:564", "Glu->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 9);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1135", "Lys->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1225", "Trp->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1181", "Gln->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:636", "Arg->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1113", "His->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1242", "Tyr->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1248", "Tyr->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to T
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1245", "Tyr->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:683", "Tyr->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:566", "Phe->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1142", "Met->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:639", "Arg->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1233", "Trp->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:556", "Asp->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:560", "Glu->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1229", "Trp->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1131", "Lys->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1226", "Trp->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:645", "Arg->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1159", "Asn->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1177", "Gln->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 8);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 8);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:638", "Arg->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1230", "Trp->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:644", "Arg->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1121", "His->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1098", "Phe->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1228", "Trp->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:680", "Tyr->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:682", "Tyr->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:568", "Phe->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1101", "Phe->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:552", "Cys->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:567", "Phe->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1207", "Thr->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1152", "Met->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:553", "Asp->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to I
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to L
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1155", "Asn->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1192", "Arg->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1129", "Xle->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1086", "Glu->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1190", "Arg->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1138", "Lys->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1184", "Gln->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1170", "Pro->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:580", "His->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to F
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1227", "Trp->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:585", "His->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1122", "His->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1243", "Tyr->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1240", "Tyr->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1114", "His->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1238", "Tyr->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1150", "Met->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1097", "Phe->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1091", "Phe->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1241", "Tyr->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1084", "Glu->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:939", "Cys->methylaminoAla", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:614", "Met->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1137", "Lys->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:630", "Gln->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:657", "Ser->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:659", "Thr->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:610", "Met->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 5);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:675", "Trp->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1081", "Glu->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:600", "Lys->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:635", "Gln->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:642", "Arg->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to C
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1143", "Met->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1125", "Xle->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to K
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:640", "Arg->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1075", "Asp->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to T
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1087", "Glu->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1191", "Arg->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:594", "Lys->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:616", "Asn->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to T
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1185", "Gln->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1078", "Glu->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:624", "Pro->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:681", "Tyr->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1132", "Lys->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:837", "Arg->Npo", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1178", "Gln->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // W to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("W");
        modificationName = "W to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1236", "Trp->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to H
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:641", "Arg->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1099", "Phe->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1095", "Phe->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1092", "Phe->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1072", "Asp->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1149", "Met->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1161", "Asn->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1096", "Phe->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1144", "Met->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to A
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:648", "Ser->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to F
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:678", "Tyr->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:548", "Cys->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1069", "Asp->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to I
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to L
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1083", "Glu->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:622", "Asn->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:595", "Lys->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to G
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to G";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:544", "Ala->Gly", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to D
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:562", "Glu->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to V
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to N
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1183", "Gln->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to S
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:658", "Thr->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1076", "Asp->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 7);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1133", "Lys->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1179", "Gln->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:617", "Asn->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:905", "Leu->MetOx", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:601", "Xle->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1067", "Asp->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1156", "Asn->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:623", "Pro->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1094", "Phe->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1195", "Arg->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:582", "His->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1119", "His->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1115", "His->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1244", "Tyr->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1120", "His->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1063", "Cys->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:662", "Thr->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1059", "Cys->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1151", "Met->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:613", "Met->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to P
        atomChainAdded = new AtomChain();
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:604", "Xle->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1145", "Met->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:664", "Thr->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1065", "Cys->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:555", "Asp->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:561", "Glu->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:563", "Glu->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:597", "Lys->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:631", "Gln->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:596", "Lys->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:621", "Asn->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:632", "Gln->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 6);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1203", "Thr->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:588", "Xle->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1082", "Glu->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:629", "Pro->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:598", "Lys->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1182", "Gln->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1126", "Xle->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:627", "Pro->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1166", "Pro->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:899", "Met->Hpg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1194", "Arg->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1080", "Glu->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 5);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1136", "Lys->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:633", "Gln->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1100", "Phe->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1116", "His->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:652", "Ser->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1062", "Cys->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1056", "Cys->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:656", "Ser->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:660", "Thr->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1073", "Asp->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1070", "Asp->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1204", "Thr->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:558", "Asp->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to A
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to A";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:571", "Gly->Ala", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1162", "Asn->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:650", "Ser->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:618", "Asn->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:589", "Xle->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1157", "Asn->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1127", "Xle->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:654", "Ser->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:540", "Ala->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:569", "Phe->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1071", "Asp->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1146", "Met->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1173", "Pro->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1160", "Asn->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.C, 0));
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1922", "Pro->HAVA", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1079", "Glu->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1134", "Lys->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1180", "Gln->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:584", "His->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:554", "Asp->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:620", "Asn->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Y to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Y");
        modificationName = "Y to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1247", "Tyr->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 4);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1064", "Cys->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:611", "Met->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1060", "Cys->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:581", "His->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:543", "Ala->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1057", "Cys->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:651", "Ser->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1209", "Thr->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:661", "Thr->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1085", "Glu->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1196", "Ser->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1205", "Thr->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:599", "Lys->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1047", "Ala->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1033", "Cys->SecNEM", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:634", "Gln->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:607", "Xle->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:590", "Xle->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1128", "Xle->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // R to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("R");
        modificationName = "R to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:637", "Arg->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:663", "Thr->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:541", "Ala->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to S
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to S";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:572", "Gly->Ser", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:626", "Pro->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1171", "Pro->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1044", "Ala->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:608", "Xle->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1168", "Pro->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1154", "Met->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1068", "Asp->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1158", "Asn->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1172", "Pro->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1089", "Glu->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1058", "Cys->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1141", "Lys->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1188", "Gln->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1208", "Thr->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:606", "Xle->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // F to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("F");
        modificationName = "F to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1102", "Phe->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:625", "Pro->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to P
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to P";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1109", "Gly->Pro", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1201", "Ser->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1199", "Ser->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1074", "Asp->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1197", "Ser->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to V
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to V";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1105", "Gly->Xle", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1163", "Asn->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1050", "Ala->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:542", "Ala->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1200", "Ser->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to T
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to T";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1111", "Gly->Thr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:547", "Cys->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to C
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to C";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:577", "Gly->Cys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1206", "Thr->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:602", "Xle->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:557", "Asp->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // H to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("H");
        modificationName = "H to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1124", "His->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.N, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:619", "Asn->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1169", "Pro->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1198", "Ser->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:551", "Cys->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // M to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("M");
        modificationName = "M to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1153", "Met->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:665", "Thr->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to I
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 8);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to I";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to L
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 8);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to L";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1051", "Ala->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to N
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to N";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1108", "Gly->Asn", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:609", "Xle->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // E to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("E");
        modificationName = "E to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1088", "Glu->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1048", "Ala->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // K to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("K");
        modificationName = "K to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1140", "Lys->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:545", "Ala->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to D
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to D";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:576", "Gly->Asp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // Q to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("Q");
        modificationName = "Q to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1187", "Gln->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0));
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:628", "Pro->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1049", "Ala->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:647", "Ser->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:550", "Cys->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1212", "Thr->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1130", "Xle->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1176", "Pro->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 2);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1046", "Ala->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:655", "Ser->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to Q
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to Q";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1110", "Gly->Gln", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // D to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0), 2);
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("D");
        modificationName = "D to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1077", "Asp->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to K
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 9);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to K";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1106", "Gly->Lys", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to E
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to E";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:574", "Gly->Glu", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // N to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("N");
        modificationName = "N to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1165", "Asn->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // I to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("I");
        modificationName = "I to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // L to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.H, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("L");
        modificationName = "L to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to M
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.S, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to M";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1107", "Gly->Met", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1045", "Ala->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:653", "Ser->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to H
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 2);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to H";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1104", "Gly->His", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // C to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.S, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("C");
        modificationName = "C to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:549", "Cys->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // T to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("T");
        modificationName = "T to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1211", "Thr->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1052", "Ala->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // V to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0));
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("V");
        modificationName = "V to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:603", "Xle->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // P to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 3);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("P");
        modificationName = "P to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1175", "Pro->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to F
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to F";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1103", "Gly->Phe", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1054", "Ala->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // S to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        atomChainRemoved.append(new AtomImpl(Atom.O, 0));
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("S");
        modificationName = "S to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:649", "Ser->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to R
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 4);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 9);
        atomChainAdded.append(new AtomImpl(Atom.N, 0), 3);
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to R";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:578", "Gly->Arg", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to Y
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 6);
        atomChainAdded.append(new AtomImpl(Atom.O, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to Y";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1112", "Gly->Tyr", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // A to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 8);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 5);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("A");
        modificationName = "A to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_TwoPlus);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:1053", "Ala->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);

        // G to W
        atomChainAdded = new AtomChain();
        atomChainAdded.append(new AtomImpl(Atom.C, 0), 9);
        atomChainAdded.append(new AtomImpl(Atom.H, 0), 7);
        atomChainAdded.append(new AtomImpl(Atom.N, 0));
        atomChainRemoved = new AtomChain();
        aminoAcidPattern = AminoAcidPattern.getAminoAcidPatternFromString("G");
        modificationName = "G to W";
        modification = new Modification(
                ModificationType.modaa, modificationName, "sub",
                atomChainAdded, atomChainRemoved, aminoAcidPattern,
                ModificationCategory.Nucleotide_Substitution_One);
        modification.setUnimodCvTerm(new CvTerm("UNIMOD", "UNIMOD:573", "Gly->Trp", null));
        defaultMods.add(modificationName);
        modificationMap.put(modificationName, modification);
    }

    /**
     * Returns an array list of the utilities modification names mapping to the
     * given PSI-MOD accession number. The accession number should be provided
     * without the leading "MOD:" part.
     *
     * @param psiModAccession the PSI-MOD accession number to look up (provided
     * without the leading "MOD:" part)
     * @return an array list of the utilities modification names mapping to the
     * given PSI-MOD accession number
     */
    public ArrayList<String> getModificationsForPsiAccession(String psiModAccession) {

        if (psiModMap.isEmpty()) {
            createPsiModMap();
        }

        return psiModMap.get(psiModAccession);
    }

    /**
     * Create the PSI-MOD map.
     */
    private void createPsiModMap() {

        Iterator<String> modificationIterator = modificationMap.keySet().iterator();

        String modName;
        while (modificationIterator.hasNext()) {

            modName = modificationIterator.next();

            CvTerm psiCvTerm = modificationMap.get(modName).getPsiModCvTerm();

            if (psiCvTerm != null) {

                String psiModAccession = psiCvTerm.getAccession();
                psiModAccession = psiModAccession.substring(psiModAccession.indexOf(":") + 1); // remove the ontology name, i.e. "MOD:"

                ArrayList<String> utilitiesModifications = psiModMap.get(psiModAccession);

                if (utilitiesModifications == null) {
                    utilitiesModifications = new ArrayList<>();
                    psiModMap.put(psiModAccession, utilitiesModifications);
                }

                utilitiesModifications.add(modName);

            }
        }

    }

    /**
     * Returns the names of all modifications in the given categories.
     *
     * @param modCategories the modification categories
     *
     * @return the names of all modifications in the given categories
     */
    public ArrayList<String> getModifications(ModificationCategory... modCategories) {

        ArrayList<String> modificationsInCategory = new ArrayList<>();

        for (String modName : modificationMap.keySet()) {

            Modification tempMod = modificationMap.get(modName);

            for (ModificationCategory modCategory : modCategories) {
                if (tempMod.getCategory() == modCategory) {
                    modificationsInCategory.add(modName);
                }
            }

        }

        return modificationsInCategory;
    }
    
    /**
     * Returns a list containing all not fixed modifications with the same
     * mass. Warning: all modifications of the profile must be loaded in the
     * modification factory.
     *
     * @param modificationMass the mass
     * @param searchParameters the search parameters
     * @return a list of all not fixed modifications with the same mass
     */
    public ArrayList<String> getSameMassNotFixedModifications(double modificationMass, SearchParameters searchParameters) {

        ArrayList<String> allNotFixed = getExpectedVariableModifications(searchParameters);

        return allNotFixed.stream()
                .filter(modName -> getModification(modName).getMass() == modificationMass)
                .collect(Collectors.toCollection(ArrayList::new));
    }

    /**
     * Returns the expected variable modifications given the search parameters.
     *
     * @param searchParameters the search parameters
     * @return expected variable modifications given the search parameters
     */
    public ArrayList<String> getExpectedVariableModifications(SearchParameters searchParameters) {

        ModificationParameters modificationParameters = searchParameters.getModificationParameters();
        ArrayList<String> expectedModifications = modificationParameters.getAllNotFixedModifications();

        MetaMorpheusParameters metaMorpheusParameters
                = (MetaMorpheusParameters) searchParameters.getIdentificationAlgorithmParameter(
                        Advocate.metaMorpheus.getIndex());

        if (metaMorpheusParameters != null && metaMorpheusParameters.runGptm()) {
            ArrayList<String> gPtmModifications = getModifications(
                    metaMorpheusParameters.getGPtmCategories().toArray(new ModificationCategory[0])
            );

            // remove fixed modifications
            for (String fixedMod : modificationParameters.getFixedModifications()) {
                gPtmModifications.remove(fixedMod);
            }

            // add to the list of expected modifications if not allready included
            for (String tempMod : gPtmModifications) {
                if (!expectedModifications.contains(tempMod)) {
                    expectedModifications.add(tempMod);
                }
            }
        }

        return expectedModifications;
    }
}
