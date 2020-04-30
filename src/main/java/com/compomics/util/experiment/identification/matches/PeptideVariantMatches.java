package com.compomics.util.experiment.identification.matches;

import com.compomics.util.db.object.DbObject;
import com.compomics.util.experiment.biology.variants.Variant;
import java.util.HashMap;

/**
 * This class represents a set of variants found on a peptide.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynski
 */
public class PeptideVariantMatches extends DbObject {

    /**
     * Empty default constructor
     */
    public PeptideVariantMatches() {
        lengthDiff = 0;
        variantMatches = null;
    }

    /**
     * The length difference induced by the variants. +1 corresponds at a
     * protein sequence of 1 amino acid longer than the peptide sequence.
     */
    private final int lengthDiff;

    /**
     * The variant matches in a map indexed by 0 based position on the peptide.
     */
    private final HashMap<Integer, Variant> variantMatches;

    /**
     * Constructor.
     *
     * @param variantMatches the variants in a map indexed by 0 based position
     * on the peptide
     * @param lengthDiff the length difference induced by the variants
     */
    public PeptideVariantMatches(
            HashMap<Integer, Variant> variantMatches, 
            int lengthDiff
    ) {

        if (variantMatches.isEmpty()) {

            throw new IllegalArgumentException("Empty map of variant provided to constructor.");

        }

        this.variantMatches = variantMatches;
        this.lengthDiff = lengthDiff;

    }

    /**
     * Returns the length difference induced by the variants. +1 corresponds at
     * a protein sequence of 1 amino acid longer than the peptide sequence.
     *
     * @return the length difference induced by the variants
     */
    public int getLengthDiff() {

        return lengthDiff;

    }

    /**
     * Returns the map of variants indexed by 0 based position on the peptide.
     *
     * @return the map of variants
     */
    public HashMap<Integer, Variant> getVariantMatches() {

        return variantMatches;

    }
}
