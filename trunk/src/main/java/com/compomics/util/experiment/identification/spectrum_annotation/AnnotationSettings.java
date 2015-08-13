package com.compomics.util.experiment.identification.spectrum_annotation;

import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.NeutralLoss;
import com.compomics.util.experiment.biology.PTM;
import com.compomics.util.experiment.biology.PTMFactory;
import com.compomics.util.experiment.identification.identification_parameters.SearchParameters;
import com.compomics.util.experiment.identification.SpectrumIdentificationAssumption;
import com.compomics.util.preferences.SequenceMatchingPreferences;
import java.io.IOException;
import java.io.Serializable;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * This class contains the spectrum annotation preferences.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class AnnotationSettings implements Serializable {

    /**
     * If true, the automatic y-axis zoom excludes background peaks. False
     * includes all peaks in the auto zoom.
     */
    private boolean yAxisZoomExcludesBackgroundPeaks = true;
    /**
     * If true, all peaks are shown, false displays the annotated peaks, and the
     * non-annotated in the background.
     */
    private boolean showAllPeaks = false;
    /**
     * The intensity limit used when only the most intense peaks are to be
     * annotated.
     */
    private double intensityLimit = 0.75;
    /**
     * Shall automatic annotation be used.
     */
    private boolean automaticAnnotation = true;
    /**
     * The types of ions to annotate.
     */
    private HashMap<Ion.IonType, HashSet<Integer>> selectedIonsMap = new HashMap<Ion.IonType, HashSet<Integer>>(4);
    /**
     * List of neutral losses to annotate.
     */
    private ArrayList<NeutralLoss> neutralLossesList;
    /**
     * If true neutral losses will be automatically deduced from the spectrum
     * identification assumption.
     */
    private Boolean neutralLossesAuto = true;
    /**
     * Fragment ion accuracy used for peak matching.
     */
    private double fragmentIonAccuracy;
    /**
     * Indicates whether the fragment ion accuracy is in ppm.
     */
    private Boolean fragmentIonPpm = false;
    /**
     * If true, the automatic forward ion de novo tags are shown.
     */
    private boolean showForwardIonDeNovoTags = false;
    /**
     * If true, the automatic rewind ion de novo tags are shown.
     */
    private boolean showRewindIonDeNovoTags = false;
    /**
     * The charge for the fragment ions in the de novo sequencing.
     */
    private int deNovoCharge = 1;
    /**
     * If there are more than one matching peak for a given annotation setting
     * this value to true results in the most accurate peak being annotated,
     * while setting this to false annotates the most intense peak.
     */
    private Boolean highResolutionAnnotation = true;

    /**
     * Constructor.
     */
    public AnnotationSettings() {
    }

    /**
     * Returns the annotation preferences specific to a spectrum and an
     * identification assumption.
     *
     * @param spectrumKey the key of the spectrum to annotate
     * @param spectrumIdentificationAssumption the spectrum identification
     * assumption to annotate with
     * @param sequenceMatchingPreferences the sequence matching preferences for
     * peptide to protein mapping
     * @param ptmSequenceMatchingPreferences the sequence matching preferences
     * for PTM to peptide mapping
     *
     * @return the annotation preferences specific to a spectrum and an
     * identification assumption
     *
     * @throws IOException exception thrown whenever an error occurred while
     * interacting with a file while mapping potential modification sites
     * @throws InterruptedException exception thrown whenever a threading issue
     * occurred while mapping potential modification sites
     * @throws ClassNotFoundException exception thrown whenever an error
     * occurred while deserializing an object from the ProteinTree
     * @throws SQLException exception thrown whenever an error occurred while
     * interacting with the ProteinTree
     */
    public SpecificAnnotationSettings getSpecificAnnotationPreferences(String spectrumKey, SpectrumIdentificationAssumption spectrumIdentificationAssumption, 
            SequenceMatchingPreferences sequenceMatchingPreferences, SequenceMatchingPreferences ptmSequenceMatchingPreferences) 
            throws IOException, InterruptedException, ClassNotFoundException, SQLException {

        SpecificAnnotationSettings specificAnnotationPreferences = new SpecificAnnotationSettings(spectrumKey, spectrumIdentificationAssumption);
        specificAnnotationPreferences.setNeutralLossesAuto(neutralLossesAuto);
        if (neutralLossesAuto) {
            specificAnnotationPreferences.setNeutralLossesMap(SpectrumAnnotator.getDefaultLosses(spectrumIdentificationAssumption, sequenceMatchingPreferences, ptmSequenceMatchingPreferences));
        } else {
            NeutralLossesMap tempNeutralLossesMap = new NeutralLossesMap();
            for (NeutralLoss neutralLoss : getNeutralLosses()) {
                tempNeutralLossesMap.addNeutralLoss(neutralLoss, 1, 1);
            }
            specificAnnotationPreferences.setNeutralLossesMap(tempNeutralLossesMap);
        }
        ArrayList<Integer> charges = new ArrayList<Integer>(4);
        int precursorCharge = spectrumIdentificationAssumption.getIdentificationCharge().value;
        if (precursorCharge == 1) {
            charges.add(precursorCharge);
        } else {
            for (int charge = 1; charge < precursorCharge; charge++) {
                charges.add(charge);
            }
        }
        specificAnnotationPreferences.setSelectedCharges(charges);
        specificAnnotationPreferences.setSelectedIonsMap((HashMap<Ion.IonType, HashSet<Integer>>) selectedIonsMap.clone());
        specificAnnotationPreferences.setFragmentIonAccuracy(fragmentIonAccuracy);
        specificAnnotationPreferences.setFragmentIonPpm(fragmentIonPpm);
        return specificAnnotationPreferences;
    }

    /**
     * Constructor setting preferences from search parameters.
     *
     * @param searchParameters the search parameters
     */
    public AnnotationSettings(SearchParameters searchParameters) {
        setPreferencesFromSearchParameters(searchParameters);
    }

    /**
     * Updates the annotation parameters based on search parameters.
     *
     * @param searchParameters the search parameters where to take the
     * information from
     */
    public void setPreferencesFromSearchParameters(SearchParameters searchParameters) {
        clearIonTypes();
        addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, searchParameters.getIonSearched1());
        addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, searchParameters.getIonSearched2());
        addIonType(Ion.IonType.TAG_FRAGMENT_ION, searchParameters.getIonSearched1());
        addIonType(Ion.IonType.TAG_FRAGMENT_ION, searchParameters.getIonSearched2());
        addIonType(Ion.IonType.PRECURSOR_ION);
        addIonType(Ion.IonType.IMMONIUM_ION);
        addIonType(Ion.IonType.REPORTER_ION);
        setFragmentIonAccuracy(searchParameters.getFragmentIonAccuracy());
        PTMFactory ptmFactory = PTMFactory.getInstance();
        for (String modName : searchParameters.getPtmSettings().getAllModifications()) {
            PTM ptm = ptmFactory.getPTM(modName);
            for (NeutralLoss neutralLoss : ptm.getNeutralLosses()) {
                addNeutralLoss(neutralLoss);
            }
        }
    }

    /**
     * Returns whether neutral losses are considered only for amino acids of
     * interest or not.
     *
     * @return a boolean indicating whether neutral losses are considered only
     * for amino acids of interest or not
     */
    public Boolean areNeutralLossesSequenceAuto() {
        return neutralLossesAuto;
    }
    
    /**
     * Sets whether neutral losses are considered only for amino acids of
     * interest or not.
     * 
     * @param neutralLossesAuto a boolean indicating whether neutral losses are considered only
     * for amino acids of interest or not
     */
    public void setNeutralLossesSequenceAuto(Boolean neutralLossesAuto) {
        this.neutralLossesAuto = neutralLossesAuto;
    }

    /**
     * Clears the considered neutral losses.
     */
    public void clearNeutralLosses() {
        neutralLossesList.clear();
    }

    /**
     * Returns the considered neutral losses.
     *
     * @return the considered neutral losses
     */
    public ArrayList<NeutralLoss> getNeutralLosses() {
        if (neutralLossesList == null) {
            neutralLossesList = new ArrayList<NeutralLoss>(2);
        }
        return neutralLossesList;
    }

    /**
     * Adds a neutral loss.
     *
     * @param neutralLoss a new neutral loss
     */
    public void addNeutralLoss(NeutralLoss neutralLoss) {
        if (neutralLossesList == null) {
            neutralLossesList = new ArrayList<NeutralLoss>(2);
        }
        boolean alreadyInList = false;
        for (NeutralLoss tempNeutralLoss : neutralLossesList) {
            if (neutralLoss.isSameAs(tempNeutralLoss)) {
                alreadyInList = true;
                break;
            }
        }
        if (!alreadyInList) {
            neutralLossesList.add(neutralLoss);
        }
    }

    /**
     * Returns the type of ions annotated.
     *
     * @return the type of ions annotated
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
            return new HashSet<Integer>();
        } else {
            return selectedIonsMap.get(Ion.IonType.PEPTIDE_FRAGMENT_ION);
        }
    }

    /**
     * Adds a new ion type and subtype to annotate.
     *
     * @param ionType a new ion type to annotate
     * @param subType the ion sub type
     */
    public void addIonType(Ion.IonType ionType, int subType) {
        if (!selectedIonsMap.containsKey(ionType)) {
            selectedIonsMap.put(ionType, new HashSet<Integer>());
        }
        this.selectedIonsMap.get(ionType).add(subType);
    }

    /**
     * Adds a new ion type to annotate. All subtypes will be annotated.
     *
     * @param ionType a new ion type to annotate
     */
    public void addIonType(Ion.IonType ionType) {
        if (!selectedIonsMap.containsKey(ionType)) {
            selectedIonsMap.put(ionType, new HashSet<Integer>());
        }
        for (int subType : Ion.getPossibleSubtypes(ionType)) {
            this.selectedIonsMap.get(ionType).add(subType);
        }
    }

    /**
     * Clears the ion types annotated.
     */
    public void clearIonTypes() {
        selectedIonsMap.clear();
    }

    /**
     * Sets whether the annotation settings should be automatically inferred.
     *
     * @param automaticAnnotation a boolean indicating whether the annotation settings should be automatically inferred
     */
    public void setAutomaticAnnotation(boolean automaticAnnotation) {
        this.automaticAnnotation = automaticAnnotation;

        if (automaticAnnotation) {
            neutralLossesAuto = true;
        }
    }

    /**
     * Returns whether the annotation settings should be automatically inferred.
     *
     * @return a boolean indicating whether the annotation settings should be automatically inferred
     */
    public boolean isAutomaticAnnotation() {
        return automaticAnnotation;
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
     * Sets the fragment ion accuracy.
     *
     * @param fragmentIonAccuracy the fragment ion accuracy
     */
    public void setFragmentIonAccuracy(double fragmentIonAccuracy) {
        this.fragmentIonAccuracy = fragmentIonAccuracy;
    }

    /**
     * Indicates whether the fragment ion accuracy is in ppm.
     *
     * @return a boolean indicating whether the fragment ion accuracy is in ppm
     */
    public boolean isFragmentIonPpm() {
        if (fragmentIonPpm == null) { // Backward compatibility
            fragmentIonPpm = false;
        }
        return fragmentIonPpm;
    }

    /**
     * Sets whether the fragment ion accuracy is in ppm.
     *
     * @param fragmentIonPpm a boolean indicating whether the fragment ion
     * accuracy is in ppm
     */
    public void setFragmentIonPpm(boolean fragmentIonPpm) {
        this.fragmentIonPpm = fragmentIonPpm;
    }

    /**
     * Returns the intensity limit. [0.0 - 1.0], where 0.0 means that all peaks
     * are considered for annotations, while 0.75 means that only the 75% most
     * intense peaks are considered for annotations.
     *
     * @return the intensityLimit
     */
    public double getAnnotationIntensityLimit() {
        return intensityLimit;
    }

    /**
     * Sets the annotation level. [0.0 - 1.0], where 0.0 means that all peaks
     * are considered for annotations, while 0.75 means that only the 75% most
     * intense peaks are considered for annotations.
     *
     * @param intensityLimit the intensityLimit to set
     */
    public void setIntensityLimit(double intensityLimit) {
        this.intensityLimit = intensityLimit;
    }

    /**
     * If true, all peaks are shown, false displays the annotated peaks, and the
     * non-annotated in the background.
     *
     * @return true if all peaks are to be shown
     */
    public boolean showAllPeaks() {
        return showAllPeaks;
    }

    /**
     * Set if all peaks or just the annotated ones are to be shown.
     *
     * @param showAllPeaks if all peaks are to be shown
     */
    public void setShowAllPeaks(boolean showAllPeaks) {
        this.showAllPeaks = showAllPeaks;
    }

    /**
     * Returns true if the automatic y-axis zoom excludes background peaks.
     * False if includes all peaks.
     *
     * @return true if the automatic y-axis zoom excludes background peaks
     */
    public boolean yAxisZoomExcludesBackgroundPeaks() {
        return yAxisZoomExcludesBackgroundPeaks;
    }

    /**
     * Set if the automatic y-axis zoom only considers the annotated peaks.
     *
     * @param yAxisZoomExcludesBackgroundPeaks if the automatic y-axis zoom only
     * considers the annotated peaks
     */
    public void setYAxisZoomExcludesBackgroundPeaks(boolean yAxisZoomExcludesBackgroundPeaks) {
        this.yAxisZoomExcludesBackgroundPeaks = yAxisZoomExcludesBackgroundPeaks;
    }

    /**
     * Returns true if the automatic forward ion de novo tags are to be shown.
     *
     * @return the showForwardIonDeNovoTags
     */
    public boolean showForwardIonDeNovoTags() {
        return showForwardIonDeNovoTags;
    }

    /**
     * Set if the automatic forward ion de novo tags are to be shown.
     *
     * @param showForwardIonDeNovoTags the showForwardIonDeNovoTags to set
     */
    public void setShowForwardIonDeNovoTags(boolean showForwardIonDeNovoTags) {
        this.showForwardIonDeNovoTags = showForwardIonDeNovoTags;
    }

    /**
     * Returns true if the automatic rewind ion de novo tags are to be shown.
     *
     * @return the showRewindIonDeNovoTags
     */
    public boolean showRewindIonDeNovoTags() {
        return showRewindIonDeNovoTags;
    }

    /**
     * Set if the automatic rewind ion de novo tags are to be shown.
     *
     * @param showRewindIonDeNovoTags the showRewindIonDeNovoTags to set
     */
    public void setShowRewindIonDeNovoTags(boolean showRewindIonDeNovoTags) {
        this.showRewindIonDeNovoTags = showRewindIonDeNovoTags;
    }

    /**
     * Returns the charge to use for the fragment ions in the de novo
     * sequencing.
     *
     * @return the charge to use for the fragment ions in the de novo sequencing
     */
    public int getDeNovoCharge() {
        return deNovoCharge;
    }

    /**
     * Set the charge to use for the fragment ions in the de novo sequencing
     *
     * @param deNovoCharge the charge to use for the fragment ions in the de
     * novo sequencing
     */
    public void setDeNovoCharge(int deNovoCharge) {
        this.deNovoCharge = deNovoCharge;
    }

    /**
     * Returns true if the peak annotation should be based on the most accurate
     * mz value, false bases the annotation on the most intense peak.
     *
     * @return the highResolutionAnnotation
     */
    public boolean isHighResolutionAnnotation() {
        if (highResolutionAnnotation == null) {
            highResolutionAnnotation = true;
        }
        return highResolutionAnnotation;
    }

    /**
     * Set if the peak annotation should be based on the most accurate mz value,
     * or on the most intense peak.
     *
     * @param highResolutionAnnotation the highResolutionAnnotation to set
     */
    public void setHighResolutionAnnotation(boolean highResolutionAnnotation) {
        this.highResolutionAnnotation = highResolutionAnnotation;
    }
    
    /**
     * Clones the settings.
     * 
     * @return a clone of this object
     */
    public AnnotationSettings clone() {
        AnnotationSettings annotationSettings = new AnnotationSettings();
        annotationSettings.setYAxisZoomExcludesBackgroundPeaks(yAxisZoomExcludesBackgroundPeaks);
        annotationSettings.setShowAllPeaks(showAllPeaks);
        annotationSettings.setIntensityLimit(intensityLimit);
        annotationSettings.setAutomaticAnnotation(automaticAnnotation);
        for (Ion.IonType ionType : selectedIonsMap.keySet()) {
            for (Integer subType : selectedIonsMap.get(ionType)) {
                annotationSettings.addIonType(ionType, subType);
            }
        }
        for (NeutralLoss neutralLoss : neutralLossesList) {
            annotationSettings.addNeutralLoss(neutralLoss);
        }
        annotationSettings.setNeutralLossesSequenceAuto(neutralLossesAuto);
        annotationSettings.setFragmentIonAccuracy(fragmentIonAccuracy);
        annotationSettings.setFragmentIonPpm(fragmentIonPpm);
        annotationSettings.setShowForwardIonDeNovoTags(showForwardIonDeNovoTags);
        annotationSettings.setShowRewindIonDeNovoTags(showRewindIonDeNovoTags);
        annotationSettings.setDeNovoCharge(deNovoCharge);
        annotationSettings.setHighResolutionAnnotation(highResolutionAnnotation);
        return annotationSettings;
    }
}
