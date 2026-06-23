package com.compomics.util.experiment.identification.matches_iterators;

import com.compomics.util.experiment.identification.matches.SpectrumMatch;

/**
 * Legacy facade for iterating over spectrum matches.
 */
public class PsmIterator {

    /**
     * The wrapped spectrum matches iterator.
     */
    private final SpectrumMatchesIterator spectrumMatchesIterator;

    /**
     * Constructor.
     *
     * @param spectrumMatchesIterator the spectrum matches iterator
     */
    public PsmIterator(SpectrumMatchesIterator spectrumMatchesIterator) {
        this.spectrumMatchesIterator = spectrumMatchesIterator;
    }

    /**
     * Returns the next spectrum match.
     *
     * @return the next spectrum match
     */
    public SpectrumMatch next() {
        return spectrumMatchesIterator.next();
    }

    /**
     * Returns the next spectrum match.
     *
     * @return the next spectrum match
     */
    public SpectrumMatch nextObject() {
        return next();
    }
}
