package com.compomics.util.test.experiment.sequences.matching;

import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidSequence;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import org.junit.Assert;
import junit.framework.TestCase;

/**
 * This class tests the matching of amino acids.
 *
 * @author Marc Vaudel
 */
public class SequenceMatchingTest extends TestCase {

    /**
     * Tests the import and the mapping of a few peptide sequences.
     *
     * @throws FileNotFoundException
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws SQLException
     * @throws InterruptedException
     */
    public void testSequenceMatchingPreferences() throws FileNotFoundException, IOException, ClassNotFoundException, SQLException, InterruptedException {

        String ref = "TESTIKTEST";
        String testAminoAcid = "TESTJKTEST";
        String testIndistinguishible = "TESTLKTEST";
        String test1Mutation = "TESTLKTETT";
        String test2Mutations = "TESTLKTETS";
        AminoAcidSequence aminoAcidSequence = new AminoAcidSequence("TESTIKTEST");

        SequenceMatchingParameters sequenceMatchingPreferences = new SequenceMatchingParameters();
        sequenceMatchingPreferences.setSequenceMatchingType(SequenceMatchingParameters.MatchingType.string);
        Assert.assertTrue(aminoAcidSequence.matches(ref, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(testAminoAcid, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(testIndistinguishible, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(test1Mutation, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(test2Mutations, sequenceMatchingPreferences));

        sequenceMatchingPreferences = new SequenceMatchingParameters();
        sequenceMatchingPreferences.setSequenceMatchingType(SequenceMatchingParameters.MatchingType.aminoAcid);
        Assert.assertTrue(aminoAcidSequence.matches(ref, sequenceMatchingPreferences));
        Assert.assertTrue(aminoAcidSequence.matches(testAminoAcid, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(testIndistinguishible, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(test1Mutation, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(test2Mutations, sequenceMatchingPreferences));

        sequenceMatchingPreferences = new SequenceMatchingParameters();
        sequenceMatchingPreferences.setSequenceMatchingType(SequenceMatchingParameters.MatchingType.indistiguishableAminoAcids);
        Assert.assertTrue(aminoAcidSequence.matches(ref, sequenceMatchingPreferences));
        Assert.assertTrue(aminoAcidSequence.matches(testAminoAcid, sequenceMatchingPreferences));
        Assert.assertTrue(aminoAcidSequence.matches(testIndistinguishible, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(test1Mutation, sequenceMatchingPreferences));
        Assert.assertTrue(!aminoAcidSequence.matches(test2Mutations, sequenceMatchingPreferences));
    }
}
