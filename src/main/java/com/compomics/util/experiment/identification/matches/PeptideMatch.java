package com.compomics.util.experiment.identification.matches;

import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.IdentificationMatch;
import java.util.Arrays;

/**
 * This class models a peptide match.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynski
 */
public class PeptideMatch extends IdentificationMatch {

    /**
     * The peptide.
     */
    private Peptide peptide;
    /**
     * The key of the match.
     */
    private long key;
    /**
     * The keys of the spectrum matches linking to this peptide match.
     */
    private long[] spectrumMatchesKeys;

    @Override
    public long getKey() {
        
        
        
        
        
        return key;
    }

    /**
     * Sets a new key for the match.
     *
     * @param newKey a new key for the match
     */
    public void setKey(long newKey) {
        
        
        
        
        
        this.key = newKey;
    }
    
    
    /**
     * Default Constructor for the peptide match.
     */
    
    
    public PeptideMatch(){
        
    }

    /**
     * Constructor for the peptide match.
     *
     * @param peptide the matching peptide
     * @param matchKey the key of the match as referenced in the identification
     * @param spectrumMatchKey the key of a spectrum match linked to this peptide
     */
    public PeptideMatch(Peptide peptide, long matchKey, long spectrumMatchKey) {
        
        
        
        
        
        this.peptide = peptide;
        this.key = matchKey;
        
        spectrumMatchesKeys = new long[1];
        spectrumMatchesKeys[0] = spectrumMatchKey;
        
    }

    /**
     * Getter for the peptide.
     *
     * @return the peptide
     */
    public Peptide getPeptide() {
        
        
        
        
        
        return peptide;
    }

    /**
     * Setter for the peptide.
     *
     * @param peptide a peptide
     */
    public void setPeptide(Peptide peptide) {
        
        
        
        
        
        this.peptide = peptide;
    }

    /**
     * Returns the keys of all spectra matched.
     *
     * @return the keys of all spectrum matches
     */
    public long[] getSpectrumMatchesKeys() {
        
        
        
        
        
        return spectrumMatchesKeys;
        
    }

    /**
     * Sets the spectrum matches keys.
     *
     * @param spectrumMatchesKeys the keys
     */
    public void setSpectrumMatchesKeys(long[] spectrumMatchesKeys) {
        
        
        
        
        
        this.spectrumMatchesKeys = spectrumMatchesKeys;
        
    }

    /**
     * Add a spectrum match key.
     *
     * @param spectrumMatchKey the key of a spectrum match
     */
    public void addSpectrumMatchKey(long spectrumMatchKey) {
        
        
        
        
        
        spectrumMatchesKeys =  Arrays.copyOf(spectrumMatchesKeys, spectrumMatchesKeys.length + 1);
        
        spectrumMatchesKeys[spectrumMatchesKeys.length - 1] = spectrumMatchKey;
        
    }

    /**
     * Returns the number of spectra matched.
     *
     * @return spectrum count
     */
    public int getSpectrumCount() {
        
        
        
        
        
        return spectrumMatchesKeys.length;
    }

    @Override
    public MatchType getType() {
        return MatchType.Peptide;
    }
}
