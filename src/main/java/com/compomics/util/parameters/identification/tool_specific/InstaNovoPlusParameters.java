package com.compomics.util.parameters.identification.tool_specific;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;

/**
 * InstaNovo+ specific parameters.
 *
 * @author CompOmics
 */
public class InstaNovoPlusParameters extends InstaNovoParameters {

    /**
     * Version number for deserialization.
     */
    static final long serialVersionUID = -7586968643672811482L;

    @Override
    public Advocate getAlgorithm() {
        return Advocate.instanovoPlus;
    }

    @Override
    public boolean equals(IdentificationAlgorithmParameter identificationAlgorithmParameter) {
        return identificationAlgorithmParameter instanceof InstaNovoPlusParameters
                && super.equals(identificationAlgorithmParameter);
    }
}
