package com.compomics.util.experiment.biology.aminoacids;

import com.compomics.util.experiment.biology.AminoAcid;
import com.compomics.util.experiment.biology.Atom;
import com.compomics.util.experiment.biology.AtomChain;
import com.compomics.util.experiment.biology.AtomImpl;

/**
 * Arginine.
 *
 * @author Marc Vaudel
 */
public class Arginine extends AminoAcid {

    /**
     * Serial number for backward compatibility.
     */
        static final long serialVersionUID = -5308475190007072857L;

    /**
     * Constructor.
     */
    public Arginine() {
        singleLetterCode = "R";
        threeLetterCode = "Arg";
        name = "Arginine";
        averageMass = 156.1857;
        monoisotopicMass = 156.101111;
        monoisotopicAtomChain = new AtomChain();
        monoisotopicAtomChain.append(new AtomImpl(Atom.C, 0), 6);
        monoisotopicAtomChain.append(new AtomImpl(Atom.H, 0), 12);
        monoisotopicAtomChain.append(new AtomImpl(Atom.N, 0), 4);
        monoisotopicAtomChain.append(new AtomImpl(Atom.O, 0), 1);
        subAminoAcidsWithoutCombination = new char[]{'R'};
        subAminoAcidsWithCombination = subAminoAcidsWithoutCombination;
        aminoAcidCombinations = new char[]{'X'};
        standardGeneticCode = new String[] {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"};
    }

    @Override
    public boolean iscombination() {
        return false;
    }
}
