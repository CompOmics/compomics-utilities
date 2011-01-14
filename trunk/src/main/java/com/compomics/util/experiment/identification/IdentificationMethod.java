package com.compomics.util.experiment.identification;

import com.compomics.util.experiment.personalization.ExperimentObject;

/**
 * This class will contain all methods used to obtain identifications.
 * User: Marc
 * Date: Nov 11, 2010
 * Time: 3:56:49 PM
 */
public class IdentificationMethod extends ExperimentObject {

    /**
     * index for identification method based on peptide mass fingerprinting
     */
    public final static int PEPTIDE_FINGERPRINTING = 0;

    /**
     * index for identification method based on MS2 fragment ion matching
     */
    public final static int MS2_IDENTIFICATION = 1;

    /**
     * index of the method
     */
    private int index;

    /**
     * Constructor for the identification method
     * @param index the index of the method as indexed by the static fields
     */
    public IdentificationMethod(int index) {
        this.index = index;
    }

    /**
     * returns the index of the identification method
     * @return the index of the identification method
     */
    public int getIndex() {
        return index;
    }
}
