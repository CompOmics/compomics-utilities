package com.compomics.util.test.parameters.identification.tool_specific;

import com.compomics.util.parameters.identification.tool_specific.InstaNovoParameters;
import junit.framework.TestCase;
import org.junit.Assert;

/**
 * Tests for InstaNovo specific parameters.
 *
 * @author CompOmics
 */
public class TestInstaNovoParameters extends TestCase {

    /**
     * Tests the desktop-oriented default batch size.
     */
    public void testDefaultBatchSize() {

        InstaNovoParameters parameters = new InstaNovoParameters();

        Assert.assertEquals(InstaNovoParameters.DEFAULT_BATCH_SIZE, parameters.getBatchSize());
        Assert.assertTrue(parameters.toString(false).contains("BATCH_SIZE=" + InstaNovoParameters.DEFAULT_BATCH_SIZE));
    }
}
