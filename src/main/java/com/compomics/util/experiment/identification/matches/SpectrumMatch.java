package com.compomics.util.experiment.identification.matches;

import com.compomics.util.experiment.identification.IdentificationMatch;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.identification.spectrum_assumptions.TagAssumption;
import com.compomics.util.experiment.personalization.ExperimentObject;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Stream;

/**
 * This class models a spectrum match.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynski
 */
public class SpectrumMatch extends IdentificationMatch {

    /**
     * The key of the match.
     */
    private long key;
    /**
     * The name of the file containing the spectrum.
     */
    private String spectrumFile;
    /**
     * The title of the spectrum.
     */
    private String spectrumtitle;
    /**
     * Map of the identification algorithm peptide assumptions: advocate number &gt;
     * score &gt; assumptions.
     */
    private HashMap<Integer, TreeMap<Double, ArrayList<PeptideAssumption>>> peptideAssumptionsMap = new HashMap<>(0);
    /**
     * The best peptide assumption.
     */
    private PeptideAssumption bestPeptideAssumption;
    /**
     * Map of the identification algorithm tag assumptions: advocate number &gt;
     * score &gt; assumptions.
     */
    private HashMap<Integer, TreeMap<Double, ArrayList<TagAssumption>>> tagAssumptionsMap = new HashMap<>(0);
    /**
     * The best tag assumption.
     */
    private TagAssumption bestTagAssumption;

    /**
     * Constructor for the spectrum match.
     */
    public SpectrumMatch() {
    }

    /**
     * Sets the peptide assumption map.
     *
     * @param peptideAssumptionsMap the peptide assumption map
     */
    public void setPeptideAssumptionMap(HashMap<Integer, TreeMap<Double, ArrayList<PeptideAssumption>>> peptideAssumptionsMap) {
        
        writeDBMode();
        
        this.peptideAssumptionsMap = peptideAssumptionsMap;
    }

    /**
     * Sets the tag assumption map.
     *
     * @param tagAssumptionsMap the tag assumption map
     */
    public void setTagAssumptionMap(HashMap<Integer, TreeMap<Double, ArrayList<TagAssumption>>> tagAssumptionsMap) {
        
        writeDBMode();
        
        this.tagAssumptionsMap = tagAssumptionsMap;
    }

    /**
     * Constructor for the spectrum match.
     *
     * @param spectrumFile The name of the file containing the spectrum.
     * @param spectrumTitle The title of the spectrum.
     */
    public SpectrumMatch(
            String spectrumFile,
            String spectrumTitle
    ) {
        
        this.spectrumFile = spectrumFile;
        this.spectrumtitle = spectrumTitle;
        this.key = getKey(spectrumFile, spectrumTitle);
        
    }
    
    /**
     * Returns a key to use for the spectrum match based on the file where the spectrum was found and its title.
     * 
     * @param spectrumFile The name of the file containing the spectrum.
     * @param spectrumTitle The title of the spectrum.
     * 
     * @return The key as long.
     */
    public static long getKey(
            String spectrumFile,
            String spectrumTitle
    ) {
        return ExperimentObject.asLong(String.join("", spectrumFile, spectrumTitle));
    }
    
    /**
     * Getter for the best peptide assumption.
     *
     * @return the best peptide assumption for the spectrum
     */
    public PeptideAssumption getBestPeptideAssumption() {
        
        readDBMode();
        
        return bestPeptideAssumption;
    }

    /**
     * Setter for the best peptide assumption.
     *
     * @param bestPeptideAssumption the best peptide assumption for the spectrum
     */
    public void setBestPeptideAssumption(PeptideAssumption bestPeptideAssumption) {
        
        writeDBMode();
        
        this.bestPeptideAssumption = bestPeptideAssumption;
    }

    /**
     * Getter for the best tag assumption.
     *
     * @return the best tag assumption for the spectrum
     */
    public TagAssumption getBestTagAssumption() {
        
        readDBMode();
        
        return bestTagAssumption;
    }

    /**
     * Setter for the best tag assumption.
     *
     * @param bestTagAssumption the best tag assumption for the spectrum
     */
    public void setBestTagAssumption(TagAssumption bestTagAssumption) {
        
        writeDBMode();
        
        this.bestTagAssumption = bestTagAssumption;
    }

    /**
     * Returns the name of the file where this spectrum was found.
     * 
     * @return The name of the file where this spectrum was found.
     */
    public String getSpectrumFile() {
        
        readDBMode();
        
        return spectrumFile;
    }
    
    /**
     * Sets the spectrum file name.
     * 
     * @param spectrumFile The spectrum file name.
     */
    public void setSpectrumFile(
            String spectrumFile 
    ) {
        this.spectrumFile = spectrumFile;
    }

    /**
     * Returns the title of the spectrum.
     * 
     * @return The title of the spectrum.
     */
    public String getSpectrumTitle() {
        
        readDBMode();
        
        return spectrumtitle;
    }
    
    /**
     * Sets the spectrum title.
     * 
     * @param spectrumTitle The spectrum title.
     */
    public void setSpectrumTitle(
            String spectrumTitle 
    ) {
        this.spectrumtitle = spectrumTitle;
    }

    @Override
    public long getKey() {
        
        readDBMode();
        
        return key;
    }

    /**
     * Returns all peptide assumptions for the specified search engine indexed by their
     * score. Null if none found.
     *
     * @param advocateId the desired advocate ID
     *
     * @return all assumptions
     */
    public TreeMap<Double, ArrayList<PeptideAssumption>> getAllPeptideAssumptions(int advocateId) {
        
        readDBMode();
        
        return peptideAssumptionsMap.get(advocateId);
    }

    /**
     * Returns all tag assumptions for the specified search engine indexed by their
     * score. Null if none found.
     *
     * @param advocateId the desired advocate ID
     *
     * @return all assumptions
     */
    public TreeMap<Double, ArrayList<TagAssumption>> getAllTagAssumptions(int advocateId) {
        
        readDBMode();
        
        return tagAssumptionsMap.get(advocateId);
    }

    /**
     * Returns a stream of all peptide assumptions
     *
     * @return a stream of all peptide assumptions
     */
    public Stream<PeptideAssumption> getAllPeptideAssumptions() {
        
        readDBMode();
        
        return peptideAssumptionsMap.values().stream()
                .flatMap(algorithmMap -> algorithmMap.values().stream())
                .flatMap(assumptionsList -> assumptionsList.stream());
    }

    /**
     * Returns a stream of all tag assumptions
     *
     * @return a stream of all tag assumptions
     */
    public Stream<TagAssumption> getAllTagAssumptions() {
        
        readDBMode();
        
        return tagAssumptionsMap.values().stream()
                .flatMap(algorithmMap -> algorithmMap.values().stream())
                .flatMap(assumptionsList -> assumptionsList.stream());
    }

    /**
     * Returns the peptide assumptions map: advocate id &gt; score &gt; list of
     * assumptions.
     *
     * @return the assumptions map
     */
    public HashMap<Integer, TreeMap<Double, ArrayList<PeptideAssumption>>> getPeptideAssumptionsMap() {
        
        readDBMode();
        
        return peptideAssumptionsMap;
    }

    /**
     * Returns the tag assumptions map: advocate id &gt; score &gt; list of
     * assumptions.
     *
     * @return the assumptions map
     */
    public HashMap<Integer, TreeMap<Double, ArrayList<TagAssumption>>> getTagAssumptionsMap() {
        
        readDBMode();
        
        return tagAssumptionsMap;
    }

    /**
     * Add a hit.
     *
     * @param advocateId the index of the advocate of the new hit
     * @param peptideAssumption the new identification assumption
     */
    public void addPeptideAssumption(int advocateId, PeptideAssumption peptideAssumption) {
        
        writeDBMode();
        
        TreeMap<Double, ArrayList<PeptideAssumption>> advocateMap = peptideAssumptionsMap.get(advocateId);
        
        if (advocateMap == null) {
            
            advocateMap = new TreeMap<>();
            peptideAssumptionsMap.put(advocateId, advocateMap);
            
        }
        
        Double score = peptideAssumption.getScore();
        ArrayList<PeptideAssumption> assumptionList = advocateMap.get(score);
        
        if (assumptionList == null) {
            
            assumptionList = new ArrayList<>(1);
            advocateMap.put(score, assumptionList);
            
        }
        
        assumptionList.add(peptideAssumption);
    }

    /**
     * Add a hit.
     *
     * @param advocateId the index of the advocate of the new hit
     * @param tagAssumption the new identification assumption
     */
    public void addTagAssumption(int advocateId, TagAssumption tagAssumption) {
        
        
        writeDBMode();
        
        
        TreeMap<Double, ArrayList<TagAssumption>> advocateMap = tagAssumptionsMap.get(advocateId);
        
        if (advocateMap == null) {
            
            advocateMap = new TreeMap<>();
            tagAssumptionsMap.put(advocateId, advocateMap);
            
        }
        
        double score = tagAssumption.getScore();
        ArrayList<TagAssumption> assumptionList = advocateMap.get(score);
        
        if (assumptionList == null) {
            
            assumptionList = new ArrayList<>(1);
            advocateMap.put(score, assumptionList);
            
        }
        
        assumptionList.add(tagAssumption);
    }

    @Override
    public MatchType getType() {
        
        
        readDBMode();
        
        
        return MatchType.Spectrum;
    }

    /**
     * Removes an assumption from the mapping.
     *
     * @param peptideAssumption the peptide assumption to remove
     */
    public void removePeptideAssumption(PeptideAssumption peptideAssumption) {
        
        writeDBMode();
        
        int se = peptideAssumption.getAdvocate();
        TreeMap<Double, ArrayList<PeptideAssumption>> algorithmMap = peptideAssumptionsMap.get(se);
        ArrayList<PeptideAssumption> assumptionsList = algorithmMap.get(peptideAssumption.getScore());
        assumptionsList.remove(peptideAssumption);
        
        if (assumptionsList.isEmpty()) {
            
            algorithmMap.remove(peptideAssumption.getScore());
            
        }
        if (algorithmMap.isEmpty()) {
            
            peptideAssumptionsMap.remove(se);
            
        }
    }

    /**
     * Removes an assumption from the mapping.
     *
     * @param tagAssumption the tag assumption to remove
     */
    public void removeTagAssumption(TagAssumption tagAssumption) {
        
        writeDBMode();
        
        int se = tagAssumption.getAdvocate();
        TreeMap<Double, ArrayList<TagAssumption>> algorithmMap = tagAssumptionsMap.get(se);
        ArrayList<TagAssumption> assumptionsList = algorithmMap.get(tagAssumption.getScore());
        assumptionsList.remove(tagAssumption);
        
        if (assumptionsList.isEmpty()) {
            
            algorithmMap.remove(tagAssumption.getScore());
            
        }
        
        if (algorithmMap.isEmpty()) {
            
            tagAssumptionsMap.remove(se);
            
        }
    }

    /**
     * Indicates whether the spectrum match contains a peptide assumption.
     *
     * @return a boolean indicating whether the spectrum match contains a peptide
     * assumption
     */
    public boolean hasPeptideAssumption() {
        
        readDBMode();
        
        return peptideAssumptionsMap.values().stream()
                .flatMap(algorithmMap -> algorithmMap.values().stream())
                .anyMatch(assumptionsList -> !assumptionsList.isEmpty());
    }

    /**
     * Indicates whether the spectrum match contains a tag assumption.
     *
     * @return a boolean indicating whether the spectrum match contains a tag
     * assumption
     */
    public boolean hasTagAssumption() {
        
        readDBMode();
        
        return tagAssumptionsMap.values().stream().flatMap(algorithmMap -> algorithmMap.values().stream())
                .anyMatch(assumptionsList -> !assumptionsList.isEmpty());
    }

    /**
     * Indicates whether the spectrum match contains a peptide assumption for
     * the given advocate (see the Advocate class).
     *
     * @param advocateId The index of the advocate
     *
     * @return a boolean indicating whether the spectrum match contains an
     * assumption for the given advocate
     */
    public boolean hasPeptideAssumption(int advocateId) {
        
        readDBMode();
        
        TreeMap<Double, ArrayList<PeptideAssumption>> algorithmIds = peptideAssumptionsMap.get(advocateId);

        return algorithmIds == null ? false : algorithmIds.values().stream().anyMatch(assumptions -> !assumptions.isEmpty());
    }

    /**
     * Indicates whether the spectrum match contains a tag assumption for
     * the given advocate (see the Advocate class).
     *
     * @param advocateId The index of the advocate
     *
     * @return a boolean indicating whether the spectrum match contains an
     * assumption for the given advocate
     */
    public boolean hasTagAssumption(int advocateId) {
        
        readDBMode();
        
        TreeMap<Double, ArrayList<TagAssumption>> algorithmIds = tagAssumptionsMap.get(advocateId);

        return algorithmIds == null ? false : algorithmIds.values().stream().anyMatch(assumptions -> !assumptions.isEmpty());
    }
}
