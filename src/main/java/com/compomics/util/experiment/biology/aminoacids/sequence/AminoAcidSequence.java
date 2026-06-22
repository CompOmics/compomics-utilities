package com.compomics.util.experiment.biology.aminoacids.sequence;

import com.compomics.util.experiment.biology.aminoacids.AminoAcid;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.Util;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.amino_acid_tags.TagComponent;
import com.compomics.util.experiment.identification.utils.ModificationUtils;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import java.nio.CharBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class represents a series of amino acids with associated modifications.
 *
 * @author Marc Vaudel
 */
public class AminoAcidSequence extends ExperimentObject implements TagComponent {

    /**
     * The sequence as string.
     */
    private String sequence;
    /**
     * The sequence as string builder.
     */
    private StringBuilder sequenceStringBuilder = null;
    /**
     * The variable modifications carried by the amino acid sequence at target
     * amino acids. 1 is the first amino acid.
     */
    private ModificationMatch[] variableModifications = null;
    /**
     * Convenience array for no modifications
     */
    private static final ModificationMatch[] noMod = new ModificationMatch[0];

    /**
     * Creates a blank sequence. All maps are null.
     */
    public AminoAcidSequence() {
    }

    /**
     * Constructor taking a sequence of amino acids as input.
     *
     * @param sequence a sequence of amino acids
     */
    public AminoAcidSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     * Constructor taking a sequence of amino acids as input.
     *
     * @param sequence a sequence of amino acids
     * @param variableModifications the variable modifications of this sequence
     */
    public AminoAcidSequence(String sequence, ModificationMatch[] variableModifications) {
        this.sequence = sequence;
        this.variableModifications = variableModifications;
    }

    /**
     * Creates a sequence from another sequence.
     *
     * @param sequence the other sequence
     */
    public AminoAcidSequence(AminoAcidSequence sequence) {

        this.sequence = sequence.getSequence();
        ModificationMatch[] modificationMatches = sequence.getVariableModifications();

        if (modificationMatches.length > 0) {

            variableModifications = Arrays.stream(modificationMatches)
                    .map(ModificationMatch::clone)
                    .toArray(ModificationMatch[]::new);

        }
    }

    /**
     * Returns the sequence as String.
     *
     * @return the sequence as String
     */
    public String getSequence() {
        

        setSequenceStringBuilder(false);

        return sequence == null ? "" : sequence;

    }

    /**
     * Returns the amino acid at the given index on the sequence in its single
     * letter code. 0 is the first amino acid.
     *
     * @param aa the index on the sequence
     *
     * @return the amino acid at the given index on the sequence in its single
     * letter code
     */
    public char charAt(int aa) {
        

        setSequenceStringBuilder(false);

        return sequence.charAt(aa);

    }

    /**
     * Returns the amino acid at the given index on the sequence. 0 is the first
     * amino acid.
     *
     * @param aa the index on the sequence
     *
     * @return the amino acid at the given index on the sequence
     */
    public AminoAcid getAminoAcidAt(int aa) {

        return AminoAcid.getAminoAcid(charAt(aa));

    }

    /**
     * Sets the sequence.
     *
     * @param aminoAcidSequence the sequence
     */
    public void setSequence(String aminoAcidSequence) {        

        sequenceStringBuilder = null;

        this.sequence = aminoAcidSequence;

    }

    /**
     * replaces the amino acid at the given position by the given amino acid
     * represented by its single letter code. 0 is the first amino acid.
     *
     * @param index the index where the amino acid should be set.
     * @param aa the amino acid to be set
     */
    public void setAaAtIndex(int index, char aa) {

        setSequenceStringBuilder(true);
        sequenceStringBuilder.setCharAt(index, aa);

    }

    /**
     * Loads the sequence in the string builder.
     */
    private void setSequenceStringBuilder(boolean stringbuilder) {
        if (stringbuilder && sequenceStringBuilder == null) {

            if (sequence != null) {

                sequenceStringBuilder = new StringBuilder(sequence);
                sequence = null;

            } else {

                sequenceStringBuilder = new StringBuilder(1);

            }

        } else if (sequence == null && sequenceStringBuilder != null) {

            sequence = sequenceStringBuilder.toString();

        }
    }

    /**
     * the sequence is kept in different formats internally. Calling this method
     * removes them from the cache.
     */
    public void emptyInternalCaches() {

        sequenceStringBuilder = null;

    }

    /**
     * Indicates whether the sequence is found in the given amino acid sequence.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return a boolean indicating whether the sequence is found in the given
     * amino acid sequence
     */
    public boolean matchesIn(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingPreferences) {

        return firstIndex(aminoAcidSequence, sequenceMatchingPreferences) >= 0;

    }

    /**
     * Indicates whether the sequence is found in the given amino acid sequence.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return a boolean indicating whether the sequence is found in the given
     * amino acid sequence
     */
    public boolean matchesIn(AminoAcidSequence aminoAcidSequence, SequenceMatchingParameters sequenceMatchingPreferences) {


        return matchesIn(aminoAcidSequence.getSequence(), sequenceMatchingPreferences);

    }

    /**
     * Indicates whether the sequence matches the given amino acid sequence in
     * size and according to the given matching preferences.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return a boolean indicating whether the sequence is found in the given
     * amino acid sequence
     */
    public boolean matches(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingPreferences) {

        return length() == aminoAcidSequence.length() && firstIndex(aminoAcidSequence, sequenceMatchingPreferences) >= 0;

    }

    /**
     * Indicates whether the sequence matches the given amino acid sequence in
     * size and according to the given matching preferences.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return a boolean indicating whether the sequence is found in the given
     * amino acid sequence
     */
    public boolean matches(AminoAcidSequence aminoAcidSequence, SequenceMatchingParameters sequenceMatchingPreferences) {

        return matches(aminoAcidSequence.getSequence(), sequenceMatchingPreferences);

    }

    /**
     * Returns the first index where the amino acid sequence is found in the
     * given sequence. -1 if not found. 0 is the first amino acid.
     *
     * @param aminoAcidSequence the amino acid sequence to look into
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return the first index where the amino acid sequence is found
     */
    public int firstIndex(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingPreferences) {

        return AminoAcidPattern.getAminoAcidPatternFromString(sequence).firstIndex(aminoAcidSequence, sequenceMatchingPreferences, 0);

    }

    /**
     * Returns the length of the sequence in amino acids.
     *
     * @return the length of the sequence in amino acids
     */
    public int length() {

        if (sequence != null) {

            return sequence.length();

        } else if (sequenceStringBuilder != null) {

            return sequenceStringBuilder.length();

        } else {

            return 0;

        }
    }

    /**
     * Appends another sequence at the end of this sequence.
     *
     * @param otherSequence the other sequence to append.
     */
    public void appendCTerm(AminoAcidSequence otherSequence) {

        setSequenceStringBuilder(true);

        int previousLength = length();
        sequenceStringBuilder.append(otherSequence.getSequence());

        ModificationMatch[] otherModifications = otherSequence.getVariableModifications();

        if (otherModifications.length > 0) {

            otherModifications = Arrays.stream(otherModifications)
                    .map(ModificationMatch::clone)
                    .peek(modificationMatch -> modificationMatch.setSite(modificationMatch.getSite() + previousLength))
                    .toArray(ModificationMatch[]::new);

            if (variableModifications == null) {

                variableModifications = otherModifications;

            } else {

                ModificationMatch[] mergedModifications = new ModificationMatch[variableModifications.length + otherModifications.length];
                System.arraycopy(variableModifications, 0, mergedModifications, 0, variableModifications.length);
                System.arraycopy(otherModifications, 0, mergedModifications, variableModifications.length, otherModifications.length);
                variableModifications = mergedModifications;

            }
        }
    }

    /**
     * Appends a series of unmodified amino acids to the sequence.
     *
     * @param otherSequence a series of unmodified amino acids represented by
     * their single letter code
     */
    public void appendCTerm(String otherSequence) {

        setSequenceStringBuilder(true);
        sequenceStringBuilder.append(otherSequence);

    }

    /**
     * Inserts another sequence in this sequence.
     *
     * @param offset the index where this sequence should be inserted, 0 is the
     * first amino acid.
     * @param otherSequence the other sequence to insert.
     */
    public void insert(int offset, AminoAcidSequence otherSequence) {

        setSequenceStringBuilder(true);
        sequenceStringBuilder.insert(0, otherSequence.getSequence());
        int otherSequenceLength = otherSequence.length();
        ModificationMatch[] otherModifications = otherSequence.getVariableModifications();

        if (otherModifications.length > 0) {

            otherModifications = Arrays.stream(otherModifications)
                    .map(ModificationMatch::clone)
                    .peek(modificationMatch -> modificationMatch.setSite(modificationMatch.getSite() + offset))
                    .toArray(ModificationMatch[]::new);

            if (variableModifications == null) {

                variableModifications = otherModifications;

            } else {

                for (ModificationMatch modificationMatch : variableModifications) {

                    int oldSite = modificationMatch.getSite();

                    if (oldSite > offset + 1) {

                        int newSite = oldSite + otherSequenceLength;
                        modificationMatch.setSite(newSite);

                    }

                }

                variableModifications = Arrays.copyOf(variableModifications, variableModifications.length + otherModifications.length);
                System.arraycopy(otherModifications, 0, variableModifications, variableModifications.length, otherModifications.length);

            }
        }
    }

    /**
     * Inserts another sequence in this sequence.
     *
     * @param offset the index where this sequence should be inserted, 0 is the
     * first amino acid.
     * @param otherSequence the other sequence to insert.
     */
    public void insert(int offset, String otherSequence) {

        setSequenceStringBuilder(true);
        sequenceStringBuilder.insert(offset, otherSequence);

    }

    /**
     * Appends another sequence at the beginning of this sequence keeping the
     * original order.
     *
     * @param otherSequence the other sequence to append.
     */
    public void appendNTerm(AminoAcidSequence otherSequence) {

        insert(0, otherSequence);

    }

    /**
     * Appends a series of unmodified amino acids to the beginning sequence
     * keeping the original order.
     *
     * @param otherSequence a series of unmodified amino acids represented by
     * their single letter code
     */
    public void appendNTerm(String otherSequence) {

        insert(0, otherSequence);

    }

    /**
     * Getter for the modifications carried by this sequence in a map: aa number
     * &gt; modification matches. 1 is the first amino acid.
     *
     * @return the modifications matches as found by the search engine
     */
    public ModificationMatch[] getVariableModifications() {

        return variableModifications == null ? noMod : variableModifications;

    }

    public ArrayList<Integer> getModificationIndexes() {
        ArrayList<Integer> result = new ArrayList<>();
        for (ModificationMatch modificationMatch : getVariableModifications()) {
            if (!result.contains(modificationMatch.getSite())) {
                result.add(modificationMatch.getSite());
            }
        }
        return result;
    }

    public ArrayList<ModificationMatch> getModificationsAt(int site) {
        ArrayList<ModificationMatch> result = new ArrayList<>();
        for (ModificationMatch modificationMatch : getVariableModifications()) {
            if (modificationMatch.getSite() == site) {
                result.add(modificationMatch);
            }
        }
        return result;
    }

    public static String getTaggedResidue(
            char aminoAcid,
            String modificationName,
            ModificationParameters modificationParameters,
            int modificationType,
            boolean useHtmlColorCoding,
            boolean useShortName
    ) {
        String[] variableModifications = new String[3];
        variableModifications[1] = modificationName;
        return ModificationUtils.getTaggedModifiedSequence(
                modificationParameters,
                Character.toString(aminoAcid),
                variableModifications,
                null,
                null,
                new String[3],
                useHtmlColorCoding,
                useShortName
        );
    }

    /**
     * Returns the variable modifications indexed by site. Modifications are
     * indexed by site as follows: N-term modifications are at index 0, C-term
     * at sequence length + 1, and amino acid at 1-based index on the sequence.
     *
     * @return the variable modifications indexed by site
     */
    public String[] getIndexedVariableModifications() {

        String[] result = new String[length() + 2];

        if (variableModifications != null) {

            for (ModificationMatch modificationMatch : variableModifications) {

                String modName = modificationMatch.getModification();
                int site = modificationMatch.getSite();

                if (result[site] == null) {

                    result[site] = modName;

                } else {

                    throw new IllegalArgumentException("Two modifications (" + result[site] + " and " + modName + ") found on sequence " + getSequence() + " at site " + site + ".");

                }
            }
        }

        return result;

    }

    /**
     * Adds a modification to one of the amino acid sequence.
     *
     * @param modificationMatch the modification match
     */
    public void addVariableModification(ModificationMatch modificationMatch) {

        variableModifications = variableModifications == null ? new ModificationMatch[1] : Arrays.copyOf(variableModifications, variableModifications.length + 1);

        variableModifications[variableModifications.length - 1] = modificationMatch;

    }

    /**
     * Sets the variable modifications.
     *
     * @param variableModifications the modifications
     */
    public void setVariableModifications(ModificationMatch[] variableModifications) {
        this.variableModifications = variableModifications;

    }

    /**
     * Returns the fixed modifications for this sequence based on the given
     * modification parameters. Modifications are returned as array of
     * modification names as they appear on the sequence. N-term modifications
     * are at index 0, C-term at index sequence length + 1, and other
     * modifications at amino acid index starting from 1. An error is thrown if
     * attempting to stack modifications. Protein modifications are not taken
     * into account.
     *
     * @param nTerm boolean indicating whether the sequence is located at the
     * n-term
     * @param cTerm boolean indicating whether the sequence is located at the
     * c-term
     * @param modificationParameters the modification parameters the
     * modification parameters
     * @param modificationsSequenceMatchingParameters the sequence matching
     * parameters to use for modifications
     *
     * @return the fixed modifications for this peptide
     */
    public String[] getFixedModifications(boolean nTerm, boolean cTerm, ModificationParameters modificationParameters, SequenceMatchingParameters modificationsSequenceMatchingParameters) {

        ModificationFactory modificationFactory = ModificationFactory.getInstance();

        String[] result = new String[sequence.length() + 2];

        for (String modName : modificationParameters.getFixedModifications()) {

            Modification modification = modificationFactory.getModification(modName);

            int[] sites = ModificationUtils.getPossibleModificationSites(this, nTerm, cTerm, modification, modificationsSequenceMatchingParameters);

            for (int site : sites) {

                if (result[site] == null) {

                    result[site] = modification.getName();

                } else {

                    throw new IllegalArgumentException("Attempting to put two fixed modifications (" + result[site] + " and " + modification.getName() + " on amino acid " + site + " of peptide " + getSequence() + ".");

                }
            }
        }

        return result;

    }

    /**
     * Indicates whether another sequence has a matching sequence. Modifications
     * are considered equal when of same mass. Modifications should be loaded in
     * the PTM factory.
     *
     * @param anotherSequence the other AminoAcidPattern
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return true if the other AminoAcidPattern targets the same sequence
     */
    public boolean isSameAs(AminoAcidSequence anotherSequence, SequenceMatchingParameters sequenceMatchingPreferences) {
        if (!isSameSequenceAndModificationStatusAs(anotherSequence, sequenceMatchingPreferences)) {

            return false;

        }

        ModificationFactory modificationFactory = ModificationFactory.getInstance();

        ModificationMatch[] modificationMatches1 = getVariableModifications();
        ModificationMatch[] modificationMatches2 = anotherSequence.getVariableModifications();

        Map<Double, HashSet<ModificationMatch>> mods1 = Arrays.stream(modificationMatches1).collect(Collectors.groupingBy(
                modificationMatch -> modificationFactory.getModification(modificationMatch.getModification()).getMass(), Collectors.toCollection(HashSet::new)));
        Map<Double, HashSet<ModificationMatch>> mods2 = Arrays.stream(modificationMatches2).collect(Collectors.groupingBy(
                modificationMatch -> modificationFactory.getModification(modificationMatch.getModification()).getMass(), Collectors.toCollection(HashSet::new)));

        for (Entry<Double, HashSet<ModificationMatch>> entry1 : mods1.entrySet()) {

            HashSet<Integer> sites1 = entry1.getValue().stream()
                    .map(ModificationMatch::getSite)
                    .collect(Collectors.toCollection(HashSet::new));
            HashSet<Integer> sites2 = mods2.get(entry1.getKey()).stream()
                    .map(ModificationMatch::getSite)
                    .collect(Collectors.toCollection(HashSet::new));

            if (!sites1.equals(sites2)) {

                return false;

            }
        }

        return true;

    }

    /**
     * Indicates whether another sequence targets the same sequence without
     * accounting for PTM localization. Modifications are considered equal when
     * of same mass. Modifications should be loaded in the PTM factory.
     *
     * @param anotherSequence the other sequence
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return true if the other AminoAcidPattern targets the same sequence
     */
    public boolean isSameSequenceAndModificationStatusAs(AminoAcidSequence anotherSequence, SequenceMatchingParameters sequenceMatchingPreferences) {

        if (!matches(anotherSequence, sequenceMatchingPreferences)) {
            return false;
        }

        ModificationMatch[] modificationMatches1 = getVariableModifications();
        ModificationMatch[] modificationMatches2 = anotherSequence.getVariableModifications();

        if (modificationMatches1.length != modificationMatches2.length) {
            return false;
        }

        ModificationFactory modificationFactory = ModificationFactory.getInstance();

        Map<Double, Long> masses1 = Arrays.stream(modificationMatches1).collect(Collectors.groupingBy(
                modificationMatch -> modificationFactory.getModification(modificationMatch.getModification()).getMass(), Collectors.counting()));

        Map<Double, Long> masses2 = Arrays.stream(modificationMatches2).collect(Collectors.groupingBy(
                modificationMatch -> modificationFactory.getModification(modificationMatch.getModification()).getMass(), Collectors.counting()));

        if (masses1.size() != masses2.size()) {
            return false;
        }

        return !masses1.entrySet().stream()
                .anyMatch(entry -> masses2.get(entry.getKey()) == null
                || !entry.getValue().equals(masses2.get(entry.getKey())));
    }

    /**
     * Returns an amino acid sequence which is a reversed version of the current
     * pattern.
     *
     * @return an amino acid sequence which is a reversed version of the current
     * pattern
     */
    public AminoAcidSequence reverse() {
        setSequenceStringBuilder(false);
        AminoAcidSequence newSequence = new AminoAcidSequence((new StringBuilder(sequence)).reverse().toString());

        ModificationMatch[] newModifications = Arrays.stream(variableModifications)
                .map(ModificationMatch::clone)
                .peek(modificationMatch -> modificationMatch.setSite(sequence.length() - modificationMatch.getSite()))
                .toArray(ModificationMatch[]::new);

        newSequence.setVariableModifications(variableModifications);

        return newSequence;

    }

    /**
     * Indicates whether the given sequence contains an amino acid which is in
     * fact a combination of amino acids.
     *
     * @param sequence the sequence of interest
     *
     * @return a boolean indicating whether the given sequence contains an amino
     * acid which is in fact a combination of amino acids
     */
    public static boolean hasCombination(String sequence) {

        return hasCombination(sequence.chars());

    }

    /**
     * Indicates whether the given sequence contains an amino acid which is in
     * fact a combination of amino acids.
     *
     * @param sequence the sequence of interest
     *
     * @return a boolean indicating whether the given sequence contains an amino
     * acid which is in fact a combination of amino acids
     */
    public static boolean hasCombination(IntStream sequence) {

        return sequence.anyMatch(aa -> AminoAcid.getAminoAcid((char) aa).iscombination());

    }

    /**
     * Indicates whether the given sequence contains an amino acid which is in
     * fact a combination of amino acids.
     *
     * @param sequence the sequence of interest
     *
     * @return a boolean indicating whether the given sequence contains an amino
     * acid which is in fact a combination of amino acids
     */
    public static boolean hasCombination(char[] sequence) {

        for (int i = 0; i < sequence.length; i++) {

            char aa = sequence[i];
            AminoAcid aminoAcid = AminoAcid.getAminoAcid(aa);

            if (aminoAcid.iscombination()) {

                return true;

            }
        }

        return false;
    }

    /**
     * Returns the minimal mass that an amino acid sequence can have taking into
     * account ambiguous amino acids.
     *
     * @param sequence a sequence of amino acids represented by their single
     * letter code
     *
     * @return the minimal mass the sequence can have
     */
    public static double getMinMass(char[] sequence) {

        double minMass = 0.0;

        for (char aa : sequence) {

            AminoAcid aminoAcid = AminoAcid.getAminoAcid(aa);

            if (aminoAcid.iscombination()) {

                char[] subAa = aminoAcid.getSubAminoAcids(false);
                aminoAcid = AminoAcid.getAminoAcid(subAa[0]);
                double minMassTemp = aminoAcid.getMonoisotopicMass();

                for (int i = 1; i < subAa.length; i++) {

                    aminoAcid = AminoAcid.getAminoAcid(subAa[i]);
                    double massTemp = aminoAcid.getMonoisotopicMass();

                    if (massTemp < minMassTemp) {

                        minMassTemp = massTemp;

                    }
                }

                minMass += minMassTemp;

            } else {

                minMass += aminoAcid.getMonoisotopicMass();

            }
        }

        return minMass;

    }

    /**
     * Returns a list of all combinations which can be created from a sequence
     * when expanding ambiguous amino acids like Xs.
     *
     * @param sequence the sequence of interest
     *
     * @return a list of all combinations which can be created from a sequence
     * when expanding ambiguous amino acids like Xs
     */
    public static ArrayList<StringBuilder> getCombinations(String sequence) {

        ArrayList<StringBuilder> newCombination, combination = new ArrayList<>(1);

        for (int i = 0; i < sequence.length(); i++) {

            newCombination = new ArrayList<>(1);
            char aa = sequence.charAt(i);
            AminoAcid aminoAcid = AminoAcid.getAminoAcid(aa);

            for (char newAa : aminoAcid.getSubAminoAcids(false)) {

                if (combination.isEmpty()) {

                    StringBuilder stringBuilder = new StringBuilder(sequence.length());
                    stringBuilder.append(newAa);
                    newCombination.add(stringBuilder);

                } else {

                    for (StringBuilder stringBuilder : combination) {

                        StringBuilder newStringBuilder = new StringBuilder(sequence.length());
                        newStringBuilder.append(stringBuilder);
                        newStringBuilder.append(newAa);
                        newCombination.add(newStringBuilder);

                    }
                }
            }

            combination = newCombination;

        }

        return combination;

    }

    @Override
    public String toString() {

        setSequenceStringBuilder(false);

        return sequence;

    }

    @Override
    public String asSequence() {

        setSequenceStringBuilder(false);

        return sequence;

    }

    @Override
    public double getMass() {
        setSequenceStringBuilder(false);

        double mass = CharBuffer.wrap(sequence.toCharArray()).chars()
                .mapToObj(aa -> AminoAcid.getAminoAcid((char) aa))
                .mapToDouble(AminoAcid::getMonoisotopicMass)
                .sum();

        if (variableModifications != null) {

            ModificationFactory modificationFactory = ModificationFactory.getInstance();
            mass += Arrays.stream(variableModifications)
                    .mapToDouble(modificationMatch -> modificationFactory.getModification(modificationMatch.getModification()).getMass())
                    .sum();

        }

        return mass;

    }

    @Override
    public boolean isSameAs(TagComponent anotherCompontent, SequenceMatchingParameters sequenceMatchingPreferences) {

        if (!(anotherCompontent instanceof AminoAcidSequence)) {

            return false;

        } else {

            AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) anotherCompontent;
            return isSameAs(aminoAcidSequence, sequenceMatchingPreferences);

        }
    }

    @Override
    public boolean isSameSequenceAndModificationStatusAs(TagComponent anotherCompontent, SequenceMatchingParameters sequenceMatchingPreferences) {

        if (!(anotherCompontent instanceof AminoAcidSequence)) {

            return false;

        } else {

            AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) anotherCompontent;
            return isSameSequenceAndModificationStatusAs(aminoAcidSequence, sequenceMatchingPreferences);

        }
    }
}
