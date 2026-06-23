package com.compomics.util.parameters.identification.search;

import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This class stores the information about the modification parameters (usage,
 * colors, names).
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class ModificationParameters extends ExperimentObject {

    /**
     * Serial version number for serialization compatibility.
     */
    static final long serialVersionUID = 342611308111304721L;
    /**
     * List of the expected fixed modifications.
     */
    private ArrayList<String> fixedModifications = new ArrayList<>(0);
    /**
     * List of the expected variable modifications.
     */
    private ArrayList<String> variableModifications = new ArrayList<>(0);
    /**
     * List of variable modifications searched during the second pass search.
     */
    private ArrayList<String> refinementVariableModifications = new ArrayList<>(0);
    /**
     * List of variable modifications searched during the second pass search.
     */
    private ArrayList<String> refinementFixedModifications = new ArrayList<>(0);
    /**
     * Mapping of the expected modification names to the color used.
     */
    private HashMap<String, Integer> colors = new HashMap<>(0);
    /**
     * Back-up mapping of the modifications for portability.
     */
    private HashMap<String, Modification> backUp = new HashMap<>(0);

    /**
     * Constructor.
     */
    public ModificationParameters() {
    }

    /**
     * Constructor creating a new Modification profile based on the given one.
     *
     * @param modificationParameters the modification profile
     */
    public ModificationParameters(ModificationParameters modificationParameters) {
        fixedModifications = modificationParameters.getFixedModifications();
        variableModifications = modificationParameters.getVariableModifications();
        refinementFixedModifications = modificationParameters.getRefinementFixedModifications();
        refinementVariableModifications = modificationParameters.getRefinementVariableModifications();
        colors = modificationParameters.getColors();
        backUp = modificationParameters.getBackedUpModifications();
    }

    /**
     * Returns the expected variable modification names included in this
     * profile.
     *
     * @return the expected variable modification names included in this profile
     */
    public ArrayList<String> getVariableModifications() {
        return variableModifications;
    }

    /**
     * Returns the searched fixed modifications names.
     *
     * @return the searched fixed modifications names
     */
    public ArrayList<String> getFixedModifications() {
        return fixedModifications;
    }

    /**
     * Return the refinement variable modifications used for the second pass
     * search.
     *
     * @return the refinement variable modifications
     */
    public ArrayList<String> getRefinementVariableModifications() {
        return refinementVariableModifications;
    }

    /**
     * Clears variable modifications.
     */
    public void clearVariableModifications() {
        variableModifications.clear();
    }

    /**
     * Clears fixed modifications.
     */
    public void clearFixedModifications() {
        fixedModifications.clear();
    }

    /**
     * Clears refinement modifications.
     */
    public void clearRefinementModifications() {
        refinementVariableModifications.clear();
    }

    /**
     * Return the refinement fixed modifications used for the second pass
     * search.
     *
     * @return the refinement fixed modifications
     */
    public ArrayList<String> getRefinementFixedModifications() {
        return refinementFixedModifications;
    }

    /**
     * Returns a list of all searched modifications.
     *
     * @return a list of all searched modifications
     */
    public ArrayList<String> getAllModifications() {

        ArrayList<String> result = Stream.concat(fixedModifications.stream(), variableModifications.stream())
                .distinct()
                .collect(Collectors.toCollection(ArrayList::new));

        refinementFixedModifications.stream().filter((ptmName) -> (!result.contains(ptmName))).forEachOrdered((ptmName) -> {
            result.add(ptmName);
        });

        refinementVariableModifications.stream().filter((ptmName) -> (!result.contains(ptmName))).forEachOrdered((ptmName) -> {
            result.add(ptmName);
        });

        return result;

    }

    /**
     * Returns a list of all searched modifications but the fixed ones. Note: to
     * be fixed a modification must be fixed during the first and second pass
     * searches.
     *
     * @return a list of all searched modifications but the fixed ones
     */
    public ArrayList<String> getAllNotFixedModifications() {

        ArrayList<String> result = Stream.concat(variableModifications.stream(), refinementVariableModifications.stream())
                .distinct()
                .collect(Collectors.toCollection(ArrayList::new));

        refinementVariableModifications.stream().filter((ptmName) -> (!result.contains(ptmName))).forEachOrdered((ptmName) -> {
            result.add(ptmName);
        });

        refinementFixedModifications.stream().filter((ptmName) -> (!fixedModifications.contains(ptmName) && !result.contains(ptmName))).forEachOrdered((ptmName) -> {
            result.add(ptmName);
        });

        return result;
    }

    /**
     * Adds a variable modification. The modification name is added in the
     * variable modifications names list and the modification is saved in the
     * back-up. In case a modification with the same name was already used it
     * will be silently overwritten.
     *
     * @param modification the modification to add
     */
    public void addVariableModification(Modification modification) {
        String modName = modification.getName();
        if (!variableModifications.contains(modName)) {
            variableModifications.add(modName);
        }
        Collections.sort(variableModifications);
        backUp.put(modName, modification);
    }

    /**
     * Adds a variable refinement modification. The modification name is added
     * in the refinement modifications names list and the modification is saved
     * in the back-up. In case a modification with the same name was already
     * used it will be silently overwritten.
     *
     * @param modification the modification to add
     */
    public void addRefinementVariableModification(Modification modification) {
        String modName = modification.getName();
        if (!refinementVariableModifications.contains(modName)) {
            refinementVariableModifications.add(modName);
        }
        Collections.sort(refinementVariableModifications);
        backUp.put(modName, modification);
    }

    /**
     * Adds a fixed refinement modification. The modification name is added in
     * the refinement modifications names list and the modification is saved in
     * the back-up. In case a modification with the same name was already used
     * it will be silently overwritten.
     *
     * @param modification the modification to add
     */
    public void addRefinementFixedModification(Modification modification) {
        String modName = modification.getName();
        if (!refinementFixedModifications.contains(modName)) {
            refinementFixedModifications.add(modName);
        }
        Collections.sort(refinementFixedModifications);
        backUp.put(modName, modification);
    }

    /**
     * Adds a fixed modification. The modification name is added in the fixed
     * modifications names list and the modification is saved in the back-up. In
     * case a modification with the same name was already used it will be
     * silently overwritten.
     *
     * @param modification the modification to add
     */
    public void addFixedModification(Modification modification) {
        String modName = modification.getName();
        if (!fixedModifications.contains(modName)) {
            fixedModifications.add(modName);
        }
        Collections.sort(fixedModifications);
        backUp.put(modName, modification);
    }

    /**
     * Sets a new color for the given expected modification.
     *
     * @param expectedModification the name of the expected modification
     * @param color the new color
     */
    public void setColor(String expectedModification, int color) {
        colors.put(expectedModification, color);
    }

    public void setColor(String expectedModification, Color color) {
        setColor(expectedModification, color.getRGB());
    }

    /**
     * Returns the color used to code the given modification.
     *
     * @param modification the name of the given expected modification
     *
     * @return the corresponding color
     */
    public int getColor(String modification) {

        Integer color = colors.get(modification);

        if (color == null) {

            ModificationFactory modificationFactory = ModificationFactory.getInstance();
            color = modificationFactory.getColor(modification);

            setColor(modification, color);

        }

        return color;
    }

    /**
     * Returns the modification colors as a map.
     *
     * @return the modifications colors as a map
     */
    public HashMap<String, Integer> getColors() {
        return colors;
    }

    /**
     * Returns the back-ed up modification with the given name.
     *
     * @param modName the name of the modification of interest
     * @return the corresponding modification. Null if not found.
     */
    public Modification getModification(String modName) {
        return backUp.get(modName);
    }

    /**
     * Returns the modifications backed-up as a map. modification name &gt;
     * modification.
     *
     * @return the modifications backed-up as a map
     */
    public HashMap<String, Modification> getBackedUpModifications() {
        return backUp;
    }

    /**
     * Removes a modification from the list of variable modifications.
     *
     * @param modificationName the name of the modification
     */
    public void removeVariableModification(String modificationName) {
        while (variableModifications.contains(modificationName)) {
            variableModifications.remove(modificationName);
        }
    }

    /**
     * Removes a modification from the list of fixed modifications.
     *
     * @param modificationName the name of the modification
     */
    public void removeFixedModification(String modificationName) {
        while (fixedModifications.contains(modificationName)) {
            fixedModifications.remove(modificationName);
        }
    }

    /**
     * Removes a variable modification from the list of refinement
     * modifications.
     *
     * @param modificationName the name of the modification
     */
    public void removeRefinementVariableModification(String modificationName) {
        while (refinementVariableModifications.contains(modificationName)) {
            refinementVariableModifications.remove(modificationName);
        }
    }

    /**
     * Removes a fixed modification from the list of refinement modifications.
     *
     * @param modificationName the name of the modification
     */
    public void removeRefinementFixedModification(String modificationName) {
        while (refinementFixedModifications.contains(modificationName)) {
            refinementFixedModifications.remove(modificationName);
        }
    }

    /**
     * Indicates whether the modification is contained in the profile, looking
     * into all modifications (fixed, variable and refinement)
     *
     * @param modificationName the name of the modification
     * @return a boolean indicating whether the modification is contained in the
     * mapping
     */
    public boolean contains(String modificationName) {
        return variableModifications.contains(modificationName)
                || fixedModifications.contains(modificationName)
                || (refinementVariableModifications != null && refinementVariableModifications.contains(modificationName))
                || (refinementFixedModifications != null && refinementFixedModifications.contains(modificationName));
    }

    /**
     * Returns true of the two profiles are identical.
     *
     * @param otherProfile the profile to compare against
     * @return true of the two profiles are identical
     */
    public boolean equals(ModificationParameters otherProfile) {

        if (otherProfile == null) {
            return false;
        }

        if (!this.getVariableModifications().equals(otherProfile.getVariableModifications())) {
            return false;
        }
        if (!this.getFixedModifications().equals(otherProfile.getFixedModifications())) {
            return false;
        }
        if (!this.getRefinementVariableModifications().equals(otherProfile.getRefinementVariableModifications())) {
            return false;
        }
        if (!this.getRefinementFixedModifications().equals(otherProfile.getRefinementFixedModifications())) {
            return false;
        }

        if (this.colors.size() != otherProfile.colors.size()) {
            return false;
        }

        Iterator<String> colorKeys = this.colors.keySet().iterator();

        while (colorKeys.hasNext()) {
            String tempKey = colorKeys.next();
            if (!otherProfile.colors.containsKey(tempKey)) {
                return false;
            }
            if (!this.colors.get(tempKey).equals(otherProfile.colors.get(tempKey))) {
                return false;
            }
        }

        if (this.backUp.size() != otherProfile.backUp.size()) {
            return false;
        }

        Iterator<String> backupKeys = this.backUp.keySet().iterator();

        while (backupKeys.hasNext()) {
            String tempKey = backupKeys.next();
            if (!otherProfile.backUp.containsKey(tempKey)) {
                return false;
            }

            // @TODO: a test for identical modifications (and not just name as above) should be added
//            if (!this.backUp.get(tempKey).equals(otherProfile.backUp.get(tempKey))) {
//                return false;
//            }
        }

        return true;
    }
}
