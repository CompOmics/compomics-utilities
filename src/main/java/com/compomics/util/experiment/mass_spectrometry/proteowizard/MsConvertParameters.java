package com.compomics.util.experiment.mass_spectrometry.proteowizard;

import java.util.HashMap;
import java.util.Set;

/**
 * The parameters to use when running msconvert.
 *
 * @author Marc Vaudel
 */
public class MsConvertParameters {

    /**
     * The format to convert to.
     */
    private ProteoWizardMsFormat msFormat;
    /**
     * Map of filters to use. Index of the filter - value of the argument
     */
    private final HashMap<Integer, String> filters = new HashMap<>(2);

    /**
     * Constructor.
     */
    public MsConvertParameters() {

    }

    /**
     * Returns the format to convert to.
     *
     * @return the format to convert to
     */
    public ProteoWizardMsFormat getMsFormat() {
        return msFormat;
    }

    /**
     * Sets the format to convert to.
     *
     * @param msFormat the format to convert to
     */
    public void setMsFormat(
            ProteoWizardMsFormat msFormat
    ) {
    
        this.msFormat = msFormat;

    }

    /**
     * Returns the index of the filters selected.
     *
     * @return the index of the filters selected
     */
    public Set<Integer> getFilters() {
        return filters.keySet();
    }

    /**
     * Adds a filter.
     *
     * @param msConvertFilterIndex the index of the filter according to the
     * MsConvertFilter enumerator.
     * @param value the value of the filter according to the filter
     * specifications.
     */
    public void addFilter(
            Integer msConvertFilterIndex, 
            String value
    ) {
    
        filters.put(msConvertFilterIndex, value);
    
    }

    /**
     * Returns the value set for a given filter.
     *
     * @param msConvertFilterIndex the index of the filter of interest.
     *
     * @return the value set by the user
     */
    public String getValue(
            int msConvertFilterIndex
    ) {

        return filters.get(msConvertFilterIndex);

    }

    /**
     * Returns the filters map, filter index - value.
     *
     * @return the filters map
     */
    public HashMap<Integer, String> getFiltersMap() {

        return filters;

    }
}
