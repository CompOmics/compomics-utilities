package com.compomics.util.experiment.identification.matches;

import com.compomics.util.experiment.personalization.ExperimentObject;

/**
 * This class models the match between theoretic PTM and identification results.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynski
 */
public class ModificationMatch extends ExperimentObject {

    /**
     * The version UID for Serialization/Deserialization compatibility.
     */
    static final long serialVersionUID = -6377699582945317984L;
    /**
     * Convenience array for no modifications.
     */
    public static final ModificationMatch[] NO_MOD = new ModificationMatch[0];
    /**
     * The modification name. The modification can be accessed via
     * the factory.
     */
    private String modification;
    /**
     * The location in the sequence. N-term modifications
     * are at index 0, C-term at index sequence length + 1, and other
     * modifications at amino acid index starting from 1. An error is thrown if
     * attempting to stack modifications.
     */
    private int modifiedSite;
    /**
     * A boolean indicating whether the modification is confidently localized
     * onto the sequence. Not applicable to fixed or terminal modifications.
     */
    private boolean confident = false;
    /**
     * A boolean indicating whether the modification is inferred from another
     * peptide. Not applicable to fixed or terminal modifications.
     */
    private boolean inferred = false;

    /**
     * Constructor for a modification match.
     *
     * @param modificationName the name of the modification.
     * @param modifiedSite the position of the modification in the sequence, 1
     * is the first residue, 0 and sequence length + 1 the N and C termini.
     */
    public ModificationMatch(String modificationName, int modifiedSite) {
        
        this.modification = modificationName;
        this.modifiedSite = modifiedSite;
        
    }

    /**
     * Default constructor for a modification match.
     */
    public ModificationMatch() {
    }

    /**
     * Getter for the theoretic PTM name.
     *
     * @return the theoretic PTM name
     */
    public String getModification() {
        
        return modification;
    }

    public String getTheoreticPtm() {
        return getModification();
    }

    /**
     * Sets the theoretic PTM.
     *
     * @param modName the theoretic PTM name
     */
    public void setModification(String modName) {
        
        this.modification = modName;
    }

    public void setTheoreticPtm(String modName) {
        setModification(modName);
    }

    /**
     * Getter for the modification site. N-term modifications
     * are at index 0, C-term at index sequence length + 1, and other
     * modifications at amino acid index starting from 1. An error is thrown if
     * attempting to stack modifications.
     *
     * @return the index of the modification in the sequence
     */
    public int getSite() {
        
        return modifiedSite;
    }

    public int getModificationSite() {
        return getSite();
    }

    /**
     * Setter for the modification site, 1 is the first amino acid.
     *
     * @param site the index of the modification in the sequence
     */
    public void setSite(int site) {
        
        this.modifiedSite = site;
    }

    /**
     * Returns a boolean indicating whether the modification is confidently
     * localized on the sequence.
     *
     * @return a boolean indicating whether the modification is confidently
     * localized on the sequence
     */
    public boolean getConfident() {
        
        return confident;
    }

    public boolean isConfident() {
        return getConfident();
    }

    public boolean isVariable() {
        return true;
    }

    /**
     * Sets whether the modification is confidently localized on the sequence.
     *
     * @param confident a boolean indicating whether the modification is
     * confidently localized on the sequence
     */
    public void setConfident(boolean confident) {
        
        this.confident = confident;
    }

    /**
     * Returns a boolean indicating whether the modification is inferred from
     * another peptide.
     *
     * @return a boolean indicating whether the modification is inferred from
     * another peptide
     */
    public boolean getInferred() {
        
        return inferred;
    }

    /**
     * Sets whether the modification is inferred from another peptide.
     *
     * @param inferred a boolean indicating whether the modification is inferred
     * from another peptide
     */
    public void setInferred(boolean inferred) {
        
        this.inferred = inferred;
    }

    /**
     * Indicates whether this modification match is the same of another one. The
     * match is only compared based on the theoretic PTM and the variability.
     * The localization and its confidence is not taken into account.
     *
     * @param anotherModificationMatch another modification match
     *
     * @return a boolean indicating whether both modification matches are the
     * same.
     */
    public boolean isSameAs(ModificationMatch anotherModificationMatch) {
        
        return modification.equals(anotherModificationMatch.getModification());
        
    }
    
    /**
     * Clones the modification match into a new match with the same attributes.
     * 
     * @return a new modification match with the same attributes
     */
    public ModificationMatch clone() {
        
        ModificationMatch newMatch = new ModificationMatch(modification, modifiedSite);
        newMatch.setConfident(confident);
        newMatch.setInferred(inferred);
        
        return newMatch;
    }
}
