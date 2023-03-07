package com.compomics.util.experiment.identification.spectrum_annotation;

import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import java.util.HashMap;
import java.util.HashSet;

/**
 * The spectrum annotation preferences specific to a spectrum and an
 * identification assumption.
 *
 * @author Marc Vaudel
 */
public class SpecificAnnotationParameters {

    /**
     * Empty default constructor
     */
    public SpecificAnnotationParameters() {
    }

    /**
     * The precursor charge.
     */
    private int precursorCharge;
    /**
     * The types of ions to annotate.
     */
    private HashMap<Ion.IonType, HashSet<Integer>> selectedIonsMap = new HashMap<>(2);
    /**
     * If true neutral losses will be automatically deduced from the spectrum
     * identification assumption.
     */
    private boolean neutralLossesAuto = true;
    /**
     * The neutral losses searched for.
     */
    private NeutralLossesMap neutralLossesMap = new NeutralLossesMap();
    /**
     * The fragment charge to be searched for.
     */
    private HashSet<Integer> selectedCharges = new HashSet<>(1);
    /**
     * Fragment ion accuracy used for peak matching.
     */
    private double fragmentIonAccuracy;
    /**
     * Indicates whether the fragment ion accuracy is in ppm.
     */
    private boolean fragmentIonPpm = false;

    /**
     * Sets the precursor charge.
     *
     * @param precursorCharge The precursor charge.
     */
    public void setPrecursorCharge(
            int precursorCharge
    ) {
        this.precursorCharge = precursorCharge;
    }

    /**
     * Returns the charge of the precursor.
     *
     * @return the charge of the precursor
     */
    public int getPrecursorCharge() {
        return precursorCharge;
    }

    /**
     * Returns the map of ions to annotate.
     *
     * @return the map of ions to annotate
     */
    public HashMap<Ion.IonType, HashSet<Integer>> getIonTypes() {
        return selectedIonsMap;
    }

    /**
     * Returns the type of peptide fragment ions annotated.
     *
     * @return the type of peptide fragment ions annotated
     */
    public HashSet<Integer> getFragmentIonTypes() {
        if (selectedIonsMap.get(Ion.IonType.PEPTIDE_FRAGMENT_ION) == null) {
            return new HashSet<>(0);
        } else {
            return selectedIonsMap.get(Ion.IonType.PEPTIDE_FRAGMENT_ION);
        }
    }

    /**
     * Sets the map of ions to annotate.
     *
     * @param selectedIonsMap the map of ions to annotate
     */
    public void setSelectedIonsMap(
            HashMap<Ion.IonType, HashSet<Integer>> selectedIonsMap
    ) {
        this.selectedIonsMap = selectedIonsMap;
    }

    /**
     * Clears the ion types annotated.
     */
    public void clearIonTypes() {
        selectedIonsMap.clear();
    }

    /**
     * Adds a new ion type and subtype to annotate.
     *
     * @param ionType a new ion type to annotate
     * @param subType the ion sub type
     */
    public void addIonType(
            Ion.IonType ionType,
            int subType
    ) {
        if (!selectedIonsMap.containsKey(ionType)) {
            selectedIonsMap.put(ionType, new HashSet<>(1));
        }
        this.selectedIonsMap.get(ionType).add(subType);
    }

    /**
     * Adds a new ion type to annotate. All subtypes will be annotated.
     *
     * @param ionType a new ion type to annotate
     */
    public void addIonType(
            Ion.IonType ionType
    ) {
        if (!selectedIonsMap.containsKey(ionType)) {
            selectedIonsMap.put(ionType, new HashSet<>(1));
        }
        for (int subType : Ion.getPossibleSubtypes(ionType)) {
            this.selectedIonsMap.get(ionType).add(subType);
        }
    }

    /**
     * Returns the map of neutral losses to annotate.
     *
     * @return the map of neutral losses to annotate
     */
    public NeutralLossesMap getNeutralLossesMap() {
        return neutralLossesMap;
    }

    /**
     * Sets the map of neutral losses to annotate.
     *
     * @param neutralLossesMap the map of neutral losses to annotate
     */
    public void setNeutralLossesMap(
            NeutralLossesMap neutralLossesMap
    ) {
        this.neutralLossesMap = neutralLossesMap;
    }

    /**
     * Clears the considered neutral losses.
     */
    public void clearNeutralLosses() {
        neutralLossesMap.clearNeutralLosses();
    }

    /**
     * Adds a neutral loss.
     *
     * @param neutralLoss a new neutral loss
     */
    public void addNeutralLoss(
            NeutralLoss neutralLoss
    ) {
        neutralLossesMap.addNeutralLoss(neutralLoss, 1, 1);
    }

    /**
     * Returns the charges selected for annotation.
     *
     * @return the charges selected for annotation
     */
    public HashSet<Integer> getSelectedCharges() {
        return selectedCharges;
    }

    /**
     * Sets the charges selected for annotation.
     *
     * @param selectedCharges the charges selected for annotation
     */
    public void setSelectedCharges(HashSet<Integer> selectedCharges) {
        this.selectedCharges = selectedCharges;
    }

    /**
     * Clears the selected charges.
     */
    public void clearCharges() {
        selectedCharges.clear();
    }

    /**
     * Add a charge to take into account when annotating the spectrum.
     *
     * @param selectedCharge a charge to take into account when annotating the
     * spectrum
     */
    public void addSelectedCharge(
            int selectedCharge
    ) {
        selectedCharges.add(selectedCharge);
    }

    /**
     * Returns the fragment ion accuracy.
     *
     * @return the fragment ion accuracy
     */
    public double getFragmentIonAccuracy() {
        return fragmentIonAccuracy;
    }

    /**
     * Returns the fragment ion accuracy in Da. If the tolerance is in ppm it
     * will be converted using the given reference mass.
     *
     * @param refMass the reference mass to use for the ppm to Da conversion
     *
     * @return the fragment ion accuracy
     */
    public double getFragmentIonAccuracyInDa(
            double refMass
    ) {
        if (fragmentIonPpm) {
            return fragmentIonAccuracy * refMass / 1000000;
        } else {
            return fragmentIonAccuracy;
        }
    }

    /**
     * Sets the fragment ion accuracy.
     *
     * @param fragmentIonAccuracy the fragment ion accuracy
     */
    public void setFragmentIonAccuracy(
            double fragmentIonAccuracy
    ) {
        this.fragmentIonAccuracy = fragmentIonAccuracy;
    }

    /**
     * Indicates whether the fragment ion accuracy is in ppm.
     *
     * @return a boolean indicating whether the fragment ion accuracy is in ppm
     */
    public boolean isFragmentIonPpm() {
        return fragmentIonPpm;
    }

    /**
     * Sets whether the fragment ion accuracy is in ppm.
     *
     * @param fragmentIonPpm a boolean indicating whether the fragment ion
     * accuracy is in ppm
     */
    public void setFragmentIonPpm(
            boolean fragmentIonPpm
    ) {
        this.fragmentIonPpm = fragmentIonPpm;
    }

    /**
     * Indicates whether neutral losses should be automatically selected.
     *
     * @return a boolean indicating whether neutral losses should be
     * automatically selected
     */
    public boolean isNeutralLossesAuto() {
        return neutralLossesAuto;
    }

    /**
     * Sets whether neutral losses should be automatically selected.
     *
     * @param neutralLossesAuto a boolean indicating whether neutral losses
     * should be automatically selected
     */
    public void setNeutralLossesAuto(boolean neutralLossesAuto) {
        this.neutralLossesAuto = neutralLossesAuto;
    }

    @Override
    public SpecificAnnotationParameters clone() {

        SpecificAnnotationParameters clone = new SpecificAnnotationParameters();
        clone.setPrecursorCharge(precursorCharge);
        clone.setFragmentIonAccuracy(getFragmentIonAccuracy());
        clone.setFragmentIonPpm(isFragmentIonPpm());
        clone.setNeutralLossesAuto(isNeutralLossesAuto());
        clone.setNeutralLossesMap(getNeutralLossesMap().clone());
        clone.setSelectedCharges(new HashSet<>(getSelectedCharges()));
        HashMap<Ion.IonType, HashSet<Integer>> currentIonsMap = getIonTypes(),
                newMap = new HashMap<>(currentIonsMap.size());
        for (Ion.IonType ionType : currentIonsMap.keySet()) {
            newMap.put(ionType, new HashSet<>(currentIonsMap.get(ionType)));
        }
        clone.setSelectedIonsMap(newMap);
        return clone;
    }
}
