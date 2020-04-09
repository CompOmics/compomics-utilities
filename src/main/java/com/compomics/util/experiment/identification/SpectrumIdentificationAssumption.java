package com.compomics.util.experiment.identification;

import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;
import com.compomics.util.experiment.biology.ions.impl.PrecursorIon;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.mass_spectrometry.spectra.Peak;
import com.compomics.util.experiment.personalization.ExperimentObject;
import java.util.ArrayList;

/**
 * Spectrum identification assumption made by an identification algorithm.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynski
 */
public abstract class SpectrumIdentificationAssumption extends ExperimentObject {

    /**
     * Empty default constructor
     */
    public SpectrumIdentificationAssumption() {
    }

    /**
     * The rank of the peptide assumption for the concerned spectrum.
     */
    protected int rank;
    /**
     * The advocate supporting this assumption.
     */
    protected int advocate;
    /**
     * The charge used for identification.
     */
    protected int identificationCharge;
    /**
     * The (advocate specific) score used to rank this assumption.
     */
    protected double score;
    /**
     * The identification file.
     */
    protected String identificationFile;
    /**
     * The raw score as provided by the identification algorithm.
     */
    protected double rawScore;
    /**
     * The individual amino acid scores. Null if not set.
     */
    protected ArrayList<double[]> aminoAcidScores = null;

    /**
     * Get the identification rank.
     *
     * @return the identification rank
     */
    public int getRank() {
        readDBMode();
        return rank;
    }

    /**
     * Set the rank of the PeptideAssumption.
     *
     * @param rank the rank of the PeptideAssumptio
     */
    public void setRank(
            int rank
    ) {
        writeDBMode();
        this.rank = rank;
    }

    /**
     * Get the used advocate.
     *
     * @return the advocate index
     */
    public int getAdvocate() {
        readDBMode();
        return advocate;
    }
    
    /**
     * Sets the advocate id.
     * 
     * @param advocate The advocate id.
     */
    public void setAdvocate(
            int advocate
    ){
        writeDBMode();
        this.advocate = advocate;
    }

    /**
     * Returns the score assigned by the advocate.
     *
     * @return the score
     */
    public double getScore() {
        readDBMode();
        return score;
    }

    /**
     * Sets the score. Note: if PsmScores.scoreRoundingDecimal is not null the
     * scored will be floored accordingly.
     *
     * @param score the score
     */
    public void setScore(
            double score
    ) {
        writeDBMode();
        this.score = score;
    }

    /**
     * Returns the identification file.
     *
     * @return the identification file
     */
    public String getIdentificationFile() {
        readDBMode();
        return identificationFile;
    }
    
    /**
     * Sets the identification file.
     * 
     * @param identificationFile The identification file.
     */
    public void setIdentificationFile(
            String identificationFile
    ){
        writeDBMode();
        this.identificationFile = identificationFile;
    }

    /**
     * Returns the charge used for identification.
     *
     * @return the charge used for identification
     */
    public int getIdentificationCharge() {
        readDBMode();
        return identificationCharge;
    }

    /**
     * Sets the identification charge.
     *
     * @param identificationCharge the identification charge
     */
    public void setIdentificationCharge(
            int identificationCharge
    ) {
        writeDBMode();
        this.identificationCharge = identificationCharge;
    }

    /**
     * Returns the theoretic mass of the given assumption.
     *
     * @return the theoretic mass of the given assumption
     */
    public abstract double getTheoreticMass();

    /**
     * Returns the theoretic m/z.
     *
     * @return the theoretic m/z
     */
    public double getTheoreticMz() {
        readDBMode();
        return (getTheoreticMass() + ElementaryIon.proton.getTheoreticMass() * identificationCharge) / identificationCharge;
    }

    /**
     * Returns the precursor mass error (in ppm or Da). Note that the value is
     * returns as (experimental mass - theoretical mass) and that negative
     * values thus can occur. The isotopic error can subtracted and retrieved by
     * the function getIsotopeNumber().
     *
     * @param measuredMZ the precursor m/z
     * @param ppm if true the error is returns in ppm, false returns the error
     * in Da
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     * 
     * @return the precursor mass error (in ppm or Da)
     */
    public double getDeltaMass(
            double measuredMZ, 
            boolean ppm, 
            int minIsotope, 
            int maxIsotope
    ) {
        
        readDBMode();
        
        IonMatch ionMatch = getPrecursorMatch(measuredMZ, 0.0);
        
        return ionMatch.getError(
                ppm, 
                minIsotope, 
                maxIsotope
        );
    }

    /**
     * Returns the precursor isotope number according to the number of protons.
     *
     * @param measuredMZ the measured m/z value
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     * 
     * @return the precursor isotope number according to the number of protons
     */
    public int getIsotopeNumber(
            double measuredMZ, 
            int minIsotope, 
            int maxIsotope
    ) {
        
        readDBMode();
        
        IonMatch ionMatch = getPrecursorMatch(measuredMZ, 0.0);
        
        return ionMatch.getIsotopeNumber(
                minIsotope, 
                maxIsotope
        );
        
    }

    /**
     * Returns the ion match.
     *
     * @param precursorMz The precursor measured m/z.
     * @param precursorIntensity The precursor measured intensity.
     * 
     * @return the ion match
     */
    public IonMatch getPrecursorMatch(
            double precursorMz,
            double precursorIntensity
    ) {
    
        readDBMode();
        
        return new IonMatch(
                precursorMz,
                precursorIntensity, 
                new PrecursorIon(
                        getTheoreticMass()
                ), 
                getIdentificationCharge()
        );
    }

    /**
     * Returns the raw score as provided by the identification algorithm.
     *
     * @return the raw score as provided by the identification algorithm
     */
    public double getRawScore() {
        
        readDBMode();
        
        return rawScore;
    
    }

    /**
     * Sets the raw score as provided by the identification algorithm.
     *
     * @param rawScore the raw score as provided by the identification algorithm
     */
    public void setRawScore(
            double rawScore
    ) {
    
        writeDBMode();
        
        this.rawScore = rawScore;
    
    }

    /**
     * Returns the amino acid scores. Null of not set. For Peptide objects the
     * list only contains one element, while for Tag objects the list will be
     * indexed in identically to the TagComponent list.
     *
     * @return the amino acid scores
     */
    public ArrayList<double[]> getAminoAcidScores() {
        readDBMode();
        return aminoAcidScores;
    }

    /**
     * Set the amino acid scores. For Peptide objects the list should only
     * contain one element, while for Tag objects the list should be indexed in
     * identical to the TagComponent list.
     *
     * @param aminoAcidScores the amino acid scores
     */
    public void setAminoAcidScores(
            ArrayList<double[]> aminoAcidScores
    ) {
    
        writeDBMode();
        
        this.aminoAcidScores = aminoAcidScores;
    
    }
}
