package com.compomics.util.experiment.biology.aminoacids.sequence;

import com.compomics.util.experiment.biology.aminoacids.AminoAcid;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters.MatchingType;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * An amino acid pattern is a sequence of amino acids. For example for trypsin:
 * Target R or K not followed by P. IMPORTANT: the index for the target residue
 * is by default 0.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynsk
 */
public class AminoAcidPattern extends ExperimentObject {

    /**
     * Serial number for backward compatibility.
     */
    static final long serialVersionUID = -2823716418631089876L;
    /**
     * Cache for the amino acids at target.
     */
    private HashSet<Character> aaAtTarget = null;
    /**
     * The length of the pattern, -1 if not set.
     */
    private int length = -1;
    /**
     * The list of targeted amino acids at a given index represented by their
     * single letter code. For trypsin: 0 &gt; {R, K} 1 &gt; {all but P}
     */
    private HashMap<Integer, ArrayList<Character>> residueTargeted = null;

    /**
     * Creates a blank pattern. All maps are null.
     */
    public AminoAcidPattern() {
        length = 0;
    }

    /**
     * Creates a pattern from another pattern.
     *
     * @param aminoAcidPattern the other pattern
     */
    public AminoAcidPattern(AminoAcidPattern aminoAcidPattern) {
        HashMap<Integer, ArrayList<Character>> otherTargets = aminoAcidPattern.getAaTargeted();
        if (otherTargets != null) {
            residueTargeted = new HashMap<>(otherTargets.size());
            for (int index : otherTargets.keySet()) {
                residueTargeted.put(index, (ArrayList<Character>) otherTargets.get(index).clone());
            }
        }
    }

    /**
     * Parses the amino acid pattern from the given string as created by the
     * toString() method.
     *
     * @param aminoAcidPatternAsString the amino acid pattern as created by the
     * toString() method
     *
     * @return the amino acid pattern
     */
    public static AminoAcidPattern getAminoAcidPatternFromString(String aminoAcidPatternAsString) {
        return getAminoAcidPatternFromString(aminoAcidPatternAsString, 0);
    }

    /**
     * Parses the amino acid pattern from the given string as created by the
     * toString() method.
     *
     * @param aminoAcidPatternAsString the amino acid pattern as created by the
     * toString() method
     * @param startIndex the start index of the pattern
     *
     * @return the amino acid pattern
     */
    public static AminoAcidPattern getAminoAcidPatternFromString(String aminoAcidPatternAsString, int startIndex) {

        AminoAcidPattern aminoAcidPattern = new AminoAcidPattern();

        // check if pattern contains brackets
        int cntOpenBrackets = 0, lastIndex = 0;
        while ((lastIndex = aminoAcidPatternAsString.indexOf("[", lastIndex)) != -1) {
            lastIndex += "[".length();
            cntOpenBrackets++;
        }
        int cntCloseBrackets = 0;
        lastIndex = 0;
        while ((lastIndex = aminoAcidPatternAsString.indexOf("]", lastIndex)) != -1) {
            lastIndex += "]".length();
            cntCloseBrackets++;
        }

        if (cntOpenBrackets != cntCloseBrackets) {
            throw new IllegalArgumentException("Number of opening and closing brackets unequal");
        }

        if (cntOpenBrackets == 0) {
            for (int i = 0; i < aminoAcidPatternAsString.length(); ++i) {
                ArrayList<Character> aminoAcids = new ArrayList<>(1);
                AminoAcid.getAminoAcid(aminoAcidPatternAsString.charAt(i));
                aminoAcids.add(aminoAcidPatternAsString.charAt(i));
                aminoAcidPattern.addModificationSite(i, aminoAcids);
            }
        } else {
            int pos = 0, siteIndex = -startIndex;
            while (pos < aminoAcidPatternAsString.length()) {
                if (aminoAcidPatternAsString.charAt(pos) == '[') {
                    int end = aminoAcidPatternAsString.indexOf("]", pos + 1);
                    ArrayList<Character> aminoAcids = new ArrayList<>(end - pos + 1);
                    for (int i = pos + 1; i < end; ++i) {
                        AminoAcid.getAminoAcid(aminoAcidPatternAsString.charAt(i));
                        aminoAcids.add(aminoAcidPatternAsString.charAt(i));
                    }
                    aminoAcidPattern.addModificationSite(siteIndex++, aminoAcids);
                    pos = end + 1;
                } else {
                    ArrayList<Character> aminoAcids = new ArrayList<>(1);
                    AminoAcid.getAminoAcid(aminoAcidPatternAsString.charAt(pos));
                    aminoAcids.add(aminoAcidPatternAsString.charAt(pos));
                    aminoAcidPattern.addModificationSite(siteIndex++, aminoAcids);
                    ++pos;
                }
            }
        }

        return aminoAcidPattern;
    }

    /**
     * Returns the map of targeted amino acids. Null if not set.
     *
     * @return the map of targeted amino acids
     */
    public HashMap<Integer, ArrayList<Character>> getAaTargeted() {
        return residueTargeted;
    }

    public ArrayList<Integer> getModificationIndexes() {
        return new ArrayList<>();
    }

    public ArrayList<ModificationMatch> getModificationsAt(int site) {
        return new ArrayList<>();
    }

    /**
     * Convenience constructor giving a list of targeted residues as input. For
     * instance (S, T, Y)
     *
     * @param targetResidues a list of targeted residues
     */
    public AminoAcidPattern(ArrayList<String> targetResidues) {

        ArrayList<Character> aminoAcids = targetResidues.stream()
                .map(aa -> aa.charAt(0))
                .collect(Collectors.toCollection(ArrayList::new));
        residueTargeted = new HashMap<>(1);
        residueTargeted.put(0, aminoAcids);
        length = 1;
    }

    /**
     * Swap two rows in the pattern. The first amino acid is 0.
     *
     * @param fromRow from row
     * @param toRow to row
     */
    public void swapRows(int fromRow, int toRow) {

        if (residueTargeted == null) {
            residueTargeted = new HashMap<>(1);
        }

        if (residueTargeted.size() < fromRow || fromRow < 0 || toRow < 0) {
            throw new IllegalArgumentException("Illegal row index: " + fromRow);
        }
        if (residueTargeted.size() < toRow || toRow < 0 || fromRow < 0) {
            throw new IllegalArgumentException("Illegal row index: " + toRow);
        }

        ArrayList<Character> toRowDataTarget = residueTargeted.get(toRow);
        residueTargeted.put(toRow, residueTargeted.get(fromRow));
        residueTargeted.put(fromRow, toRowDataTarget);

        // TODO: if an error should occur, an index shifting should be added here
        aaAtTarget = null;
    }

    /**
     * Returns the index of the amino acid of interest in the pattern. Null if
     * none.
     *
     * @return the index of the amino acid of interest in the pattern.
     */
    public int getTarget() {

        return 0;
    }

    /**
     * Returns the minimal index where amino acids are found.
     *
     * @return the minimal index where amino acids are found
     */
    public int getMinIndex() {

        return Math.min(residueTargeted.keySet().stream()
                .mapToInt(index -> index)
                .min()
                .orElse(0),
                 0);
    }

    /**
     * Returns the maximal index where amino acids are found.
     *
     * @return the maximal index where amino acids are found
     */
    public int getMaxIndex() {

        return Math.max(residueTargeted.keySet().stream()
                .mapToInt(index -> (int) index)
                .max()
                .orElse(0),
                 0);
    }

    /**
     * Sets the index of the amino acid of interest in the pattern.
     *
     * @param target the index of the amino acid of interest in the pattern.
     */
    public void setTarget(Integer target) {

        if (residueTargeted == null) {
            residueTargeted = new HashMap<>(1);
        }

        if (residueTargeted.size() > 0 && !residueTargeted.containsKey(target)) {
            throw new IllegalArgumentException("Target number exceeds residue site for index shifting.");
        }

        HashMap<Integer, ArrayList<Character>> residueTargetedTmp = new HashMap<>();
        for (HashMap.Entry<Integer, ArrayList<Character>> entry : residueTargeted.entrySet()) {
            residueTargetedTmp.put(entry.getKey() - target, entry.getValue());
        }
        residueTargeted = residueTargetedTmp;
        aaAtTarget = null;
    }

    /**
     * Returns the targeted amino acids at position "target". An empty list if
     * none.
     *
     * @return the targeted amino acids at position "target"
     */
    public ArrayList<Character> getAminoAcidsAtTarget() {
        return getTargetedAA(0);
    }

    /**
     * Returns a set containing the amino acids at target.
     *
     * @return a set containing the amino acids at target
     */
    public HashSet<Character> getAminoAcidsAtTargetSet() {
        if (aaAtTarget == null) {
            ArrayList<Character> aaAtTargetList = getAminoAcidsAtTarget();
            aaAtTarget = new HashSet<>(aaAtTargetList);
        }
        return aaAtTarget;
    }

    /**
     * Sets the amino acids targeted at a given index. The first amino acid is
     * 0. Previous value will be silently overwritten.
     *
     * @param index the index in the pattern
     * @param targets the amino acids targeted
     */
    public void setTargeted(int index, ArrayList<Character> targets) {
        if (residueTargeted == null) {
            residueTargeted = new HashMap<>(1);
        }
        residueTargeted.put(index, targets);
        if (index + 1 > length) {
            length = index + 1;
        }
        aaAtTarget = null;
    }

    /**
     * Excludes the given amino acids from the targeted amino acids at the given
     * index.
     *
     * @param index the index of the excluded amino acid
     * @param exceptions the amino acids to exclude
     */
    public void setExcluded(int index, ArrayList<Character> exceptions) {
        if (residueTargeted == null) {
            residueTargeted = new HashMap<>(1);
        }
        if (exceptions == null || exceptions.isEmpty()) {
            residueTargeted.put(index, new ArrayList<>());
        } else {
            ArrayList<Character> notExcluded = new ArrayList<>();
            ArrayList<Character> targeted = residueTargeted.get(index);
            if (targeted == null || targeted.isEmpty()) {
                for (char aa : AminoAcid.getUniqueAminoAcids()) {
                    if (!exceptions.contains(aa)) {
                        notExcluded.add(aa);
                    }
                }
            } else {
                for (Character aminoAcid : targeted) {
                    if (!exceptions.contains(aminoAcid)) {
                        notExcluded.add(aminoAcid);
                    }
                }
            }
            residueTargeted.put(index, notExcluded);
        }
        if (index + 1 > length) {
            length = index + 1;
        }
        aaAtTarget = null;
    }

    /**
     * Returns the targeted amino acids at a given index in the pattern. The
     * first amino acid is 0.
     *
     * @param index the index in the pattern
     * @return the targeted amino acids
     */
    public ArrayList<Character> getTargetedAA(int index) {
        if (residueTargeted != null) {
            ArrayList<Character> result = residueTargeted.get(index);
            if (result != null) {
                return result;
            }
        }
        return new ArrayList<>(0);
    }

    /**
     * Returns the number of targeted amino acids at the given index. The first
     * amino acid is 0.
     *
     * @param index the index of interest
     *
     * @return the number of excluded amino acids
     */
    public int getNTargetedAA(int index) {
        if (residueTargeted == null) {
            return 0;
        }
        ArrayList<Character> aas = getTargetedAA(index);
        return aas.size();
    }

    /**
     * Removes an amino acid index from the pattern. The first amino acid is 0.
     *
     * @param index the index of the amino acid to remove
     */
    public void removeAA(int index) { // @TODO: delete does not always work...

        if (residueTargeted != null) {
            ArrayList<Integer> indexes = new ArrayList<>(residueTargeted.keySet());
            Collections.sort(indexes);
            for (int aa : indexes) {
                if (aa >= index) {
                    if (aa > index) {
                        residueTargeted.put(aa - 1, residueTargeted.get(aa));
                    }
                    residueTargeted.remove(aa);
                }
            }
        }

        aaAtTarget = null;
        length = -1;
    }

    /**
     * Returns the amino acid pattern as case insensitive pattern for String
     * matching.
     *
     * @param sequenceMatchingParameters the sequence matching preferences
     * @param includeMutations if true mutated amino acids will be included
     *
     * @return the amino acid pattern as java string pattern
     */
    public Pattern getAsStringPattern(SequenceMatchingParameters sequenceMatchingParameters, boolean includeMutations) {

        MatchingType matchingType = sequenceMatchingParameters.getSequenceMatchingType();

        int tempLength = length();
        StringBuilder regexBuilder = new StringBuilder(tempLength);

        for (int i = 0; i < tempLength; i++) {

            ArrayList<Character> toAdd = new ArrayList<>(1);

            if (residueTargeted != null) {
                ArrayList<Character> tempTarget = residueTargeted.get(i);

                if (tempTarget == null || tempTarget.isEmpty()) {
                    toAdd.ensureCapacity(AminoAcid.getUniqueAminoAcids().length);
                    for (Character aa : AminoAcid.getUniqueAminoAcids()) {
                        toAdd.add(aa);
                    }
                } else {
                    for (Character aa : tempTarget) {
                        if (!toAdd.contains(aa)) {
                            toAdd.add(aa);
                        }
                        if (matchingType == MatchingType.aminoAcid || matchingType == MatchingType.indistiguishableAminoAcids) {
                            AminoAcid aminoAcid = AminoAcid.getAminoAcid(aa);
                            for (char tempAa : aminoAcid.getSubAminoAcids()) {
                                if (!toAdd.contains(tempAa)) {
                                    toAdd.add(tempAa);
                                }
                            }
                            for (char tempAa : aminoAcid.getCombinations()) {
                                if (!toAdd.contains(tempAa)) {
                                    toAdd.add(tempAa);
                                }
                            }
                            if (matchingType == MatchingType.indistiguishableAminoAcids
                                    && (aminoAcid == AminoAcid.I || aminoAcid == AminoAcid.J || aminoAcid == AminoAcid.L)) {
                                if (!toAdd.contains('I')) {
                                    toAdd.add('I');
                                }
                                if (!toAdd.contains('J')) {
                                    toAdd.add('L');
                                }
                            }
                        }
                    }
                }
            }

            Collections.sort(toAdd);

            regexBuilder.ensureCapacity(toAdd.size() + 2);
            regexBuilder.append("[");

            for (Character aa : toAdd) {
                regexBuilder.append(aa);
            }

            regexBuilder.append("]");
        }

        return Pattern.compile(regexBuilder.toString(), Pattern.CASE_INSENSITIVE);
    }

    /**
     * Returns the pattern in the PROSITE format.
     *
     * @return the pattern in the PROSITE format
     */
    public String getPrositeFormat() {

        StringBuilder result = new StringBuilder();
        int cpt = 0;

        for (int i = 0; i < length(); i++) {

                ArrayList<Character> targetedAas = getTargetedAA(i);

                if (targetedAas.isEmpty()) {
                    cpt++;
                } else if (targetedAas.size() > 15) {
                    ArrayList<Character> excludedAas = new ArrayList<>();
                    for (char aa : AminoAcid.getUniqueAminoAcids()) {
                        if (!targetedAas.contains(aa)) {
                            excludedAas.add(aa);
                        }
                    }
                    if (cpt > 0) {
                        result.append("(").append(cpt).append(")");
                        cpt = 0;
                    }
                    result.append("{");
                    for (Character aa : excludedAas) {
                        result.append(aa);
                    }
                    result.append("}");
                } else {
                    if (cpt > 0) {
                        result.append("(").append(cpt).append(")");
                        cpt = 0;
                    }
                    if (!targetedAas.isEmpty()) {
                        result.append("[");
                        for (Character aa : targetedAas) {
                            result.append(aa);
                        }
                        result.append("]");
                    }
                }

                if (i == 0) {
                    result.append("!");
                }
            }

        return result.toString();
    }

    /**
     * Returns the indexes where the amino acid pattern was found in the input.
     * 1 is the first amino acid.
     *
     * @param input the amino acid input sequence as string
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a list of indexes where the amino acid pattern was found
     */
    public int[] getIndexes(String input, SequenceMatchingParameters sequenceMatchingParameters) {

        ArrayList<Integer> result = new ArrayList<>(1);
        int index = 0;

        while ((index = firstIndex(input, sequenceMatchingParameters, index)) >= 0) {

            result.add(index + 1);
            index++;

        }

        return result.stream().mapToInt(a -> a).toArray();
    }

    /**
     * Returns the indexes where the amino acid pattern was found in the input.
     * 1 is the first amino acid.
     *
     * @param input the amino acid input sequence as AminoAcidPattern
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a list of indexes where the amino acid pattern was found
     */
    public ArrayList<Integer> getIndexes(AminoAcidPattern input, SequenceMatchingParameters sequenceMatchingParameters) {
        ArrayList<Integer> result = new ArrayList<>(1);
        int index = 0;
        while ((index = firstIndex(input, sequenceMatchingParameters, index)) >= 0) {
            result.add(index + 1);
            index++;
        }
        return result;
    }

    /**
     * Returns the first index where the amino acid pattern is found. -1 if not
     * found. 0 is the first amino acid.
     *
     * @param aminoAcidSequence the amino acid sequence to look into
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return the first index where the amino acid pattern is found
     */
    public int firstIndex(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        return firstIndex(aminoAcidSequence, sequenceMatchingParameters, 0);
    }

    /**
     * Returns the first index where the amino acid pattern is found. -1 if not
     * found. 0 is the first amino acid.
     *
     * @param aminoAcidSequence the amino acid sequence to look into
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return the first index where the amino acid pattern is found
     */
    public int firstIndex(AminoAcidSequence aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        return firstIndex(aminoAcidSequence.getSequence(), sequenceMatchingParameters, 0);
    }

    /**
     * Returns the first index where the amino acid pattern is found. -1 if not
     * found. 0 is the first amino acid.
     *
     * @param aminoAcidPattern the amino acid sequence to look into
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return the first index where the amino acid pattern is found
     */
    public int firstIndex(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        return firstIndex(aminoAcidPattern, sequenceMatchingParameters, 0);
    }

    /**
     * Indicates whether the pattern contains a subsequence of amino acids.
     *
     * @param aminoAcidSequence the amino acid sequence to look for
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return the first index where the amino acid pattern is found
     */
    public boolean contains(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        AminoAcidPattern pattern = getAminoAcidPatternFromString(aminoAcidSequence);
        return pattern.firstIndex(this, sequenceMatchingParameters) >= 0;
    }

    /**
     * Indicates whether the pattern contains a subsequence of amino acids.
     *
     * @param aminoAcidPattern the amino acid sequence to look for
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return the first index where the amino acid pattern is found
     */
    public boolean contains(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        return aminoAcidPattern.firstIndex(this, sequenceMatchingParameters) >= 0;
    }

    /**
     * Returns the first index where the amino acid pattern is found. -1 if not
     * found. 0 is the first amino acid.
     *
     * @param aminoAcidSequence the amino acid sequence to look into
     * @param sequenceMatchingParameters the sequence matching preferences
     * @param startIndex the start index where to start looking for
     *
     * @return the first index where the amino acid pattern is found
     */
    public int firstIndex(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters, int startIndex) {

        int patternLength = length();
        int aminoAcidPatternLength = aminoAcidSequence.length();
        int lastIndex = aminoAcidPatternLength - patternLength;

        for (int i = startIndex; i <= lastIndex; i++) {

            boolean match = true;

            for (int j = 0; j < patternLength; j++) {
                char aa = aminoAcidSequence.charAt(i + j);
                if (!isTargeted(aa, j, sequenceMatchingParameters)) {
                    match = false;
                    break;
                }
            }

            if (match) {
                return i;
            }
        }

        return -1;
    }

    /**
     * Returns the first index where the amino acid pattern is found in the
     * given pattern. -1 if not found. 0 is the first amino acid.
     *
     * @param aminoAcidPattern the amino acid sequence to look into
     * @param sequenceMatchingParameters the sequence matching preferences
     * @param startIndex the start index where to start looking for
     *
     * @return the first index where the amino acid pattern is found
     */
    public int firstIndex(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters, int startIndex) {

        int patternLength = length();
        int aminoAcidPatternLength = aminoAcidPattern.length();
        int lastIndex = aminoAcidPatternLength - patternLength;

        for (int i = startIndex; i <= lastIndex; i++) {
            boolean match = true;
            for (int j = 0; j < patternLength; j++) {
                ArrayList<Character> aminoAcids = aminoAcidPattern.getTargetedAA(i + j);
                if (!aminoAcids.isEmpty()) {
                    boolean aaMatched = false;
                    for (Character aa : aminoAcids) {
                        if (isTargeted(aa, j, sequenceMatchingParameters)) {
                            aaMatched = true;
                            break;
                        }
                    }
                    if (!aaMatched) {
                        match = false;
                        break;
                    }
                }
            }
            if (match) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Indicates whether the given amino acid at the given index of the pattern
     * is targeted without accounting for mutations.
     *
     * @param aa the amino acid as character
     * @param index the index in the pattern
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return true if the given amino acid at the given index of the pattern is
     * targeted
     */
    public boolean isTargeted(Character aa, int index, SequenceMatchingParameters sequenceMatchingParameters) {

        if (residueTargeted != null) {

            MatchingType matchingType = sequenceMatchingParameters.getSequenceMatchingType();
            ArrayList<Character> aaList = residueTargeted.get(index);

            if (aaList != null && !aaList.isEmpty()) {

                for (int i = 0; i < aaList.size(); i++) {
                    
                    Character targetedAA = aaList.get(i);
                    
                    if (aa.equals(targetedAA)) {
                        
                        return true;
                    
                    } else if (matchingType == MatchingType.aminoAcid || matchingType == MatchingType.indistiguishableAminoAcids) {

                        AminoAcid targetedAminoAcid = AminoAcid.getAminoAcid(targetedAA);

                        for (Character tempAA : targetedAminoAcid.getSubAminoAcids()) {
                            if (aa.equals(tempAA)) {
                                return true;
                            }
                        }

                        for (Character tempAA : targetedAminoAcid.getCombinations()) {
                            if (aa.equals(tempAA)) {
                                return true;
                            }
                        }

                        if (matchingType == MatchingType.indistiguishableAminoAcids
                                && (targetedAminoAcid == AminoAcid.I || targetedAminoAcid == AminoAcid.J || targetedAminoAcid == AminoAcid.L)) {
                            
                            if (aa == 'I' || aa == 'J' || aa == 'L') {
                            
                                return true;
                            
                            }
                        }
                    }
                }
            } else if (aaList != null) {
                return true;
            }
        }

        return false;
    }

    /**
     * Indicates whether the pattern is found in the given amino acid sequence.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the pattern is found in the given
     * amino acid sequence
     */
    public boolean matchesIn(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        return firstIndex(aminoAcidSequence, sequenceMatchingParameters) >= 0;
    }

    /**
     * Indicates whether the pattern is found in the given amino acid sequence.
     *
     * @param aminoAcidPattern the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the pattern is found in the given
     * amino acid sequence
     */
    public boolean matchesIn(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        return firstIndex(aminoAcidPattern, sequenceMatchingParameters) >= 0;
    }

    /**
     * Indicates whether the pattern is found in the given amino acid sequence
     * at the given index, where 0 is the first amino acid. Returns false if the
     * entire pattern cannot be mapped to the sequence.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     * @param index the index at which the matching should be done
     *
     * @return a boolean indicating whether the pattern is found in the given
     * amino acid sequence at the given index
     */
    public boolean matchesAt(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters, int index) {
        int startIndex = index;
        int endIndex = length();
        if (startIndex < 0) {
            return false;
        }
        if (endIndex >= aminoAcidSequence.length()) {
            return false;
        }
        String subSequence = aminoAcidSequence.substring(index, index + length());
        return matches(subSequence, sequenceMatchingParameters);
    }

    /**
     * Indicates whether the pattern matches the given amino acid sequence.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the pattern matches the given amino
     * acid sequence
     */
    public boolean matches(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        return length() == aminoAcidSequence.length() && firstIndex(aminoAcidSequence, sequenceMatchingParameters) >= 0;
    }

    /**
     * Indicates whether the pattern matches the given amino acid sequence
     *
     * @param aminoAcidPattern the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the pattern matches the given amino
     * acid sequence
     */
    public boolean matches(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        return length() == aminoAcidPattern.length() && firstIndex(aminoAcidPattern, sequenceMatchingParameters) >= 0;
    }

    /**
     * Indicates whether the given amino acid sequence starts with the pattern.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the given amino acid sequence starts
     * with the pattern
     */
    public boolean isStarting(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        int patternLength = length();
        return matchesIn(aminoAcidSequence.substring(0, patternLength), sequenceMatchingParameters);
    }

    /**
     * Indicates whether the given amino acid sequence starts with the pattern.
     *
     * @param aminoAcidPattern the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the given amino acid sequence starts
     * with the pattern
     */
    public boolean isStarting(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        int patternLength = length();
        return matchesIn(aminoAcidPattern.getSubPattern(0, patternLength, false), sequenceMatchingParameters);
    }

    /**
     * Indicates whether the given amino acid sequence ends with the pattern.
     *
     * @param aminoAcidPattern the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the given amino acid sequence ends
     * with the pattern
     */
    public boolean isEnding(AminoAcidPattern aminoAcidPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        int patternLength = length();
        return matchesIn(aminoAcidPattern.getSubPattern(aminoAcidPattern.length() - patternLength, false), sequenceMatchingParameters);
    }

    /**
     * Indicates whether the given amino acid sequence ends with the pattern.
     *
     * @param aminoAcidSequence the amino acid sequence
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return a boolean indicating whether the given amino acid sequence ends
     * with the pattern
     */
    public boolean isEnding(String aminoAcidSequence, SequenceMatchingParameters sequenceMatchingParameters) {
        int patternLength = length();
        return matchesIn(aminoAcidSequence.substring(aminoAcidSequence.length() - patternLength), sequenceMatchingParameters);
    }

    /**
     * Indicates whether another AminoAcidPattern targets the same pattern.
     * Modifications are considered equal when of same mass. Modifications
     * should be loaded in the Modification factory.
     *
     * @param anotherPattern the other AminoAcidPattern
     * @param sequenceMatchingParameters the sequence matching preferences
     *
     * @return true if the other AminoAcidPattern targets the same pattern
     */
    public boolean isSameAs(AminoAcidPattern anotherPattern, SequenceMatchingParameters sequenceMatchingParameters) {
        return matches(anotherPattern, sequenceMatchingParameters);
    }

    /**
     * Returns the length of the pattern in amino acids.
     *
     * @return the length of the pattern in amino acids
     */
    public int length() {
        if (residueTargeted == null) {
            return 0;
        }
        return residueTargeted.size();
    }

    /**
     * Computes a pattern which can be searched by standard search engines,
     * i.e., a pattern targeting a single amino acid and not a complex pattern.
     *
     * @return a pattern which can be searched by standard search engines
     */
    public AminoAcidPattern getStandardSearchPattern() {
        AminoAcidPattern result = new AminoAcidPattern();
        result.setTarget(0);
        result.setTargeted(0, getAminoAcidsAtTarget());
        return result;
    }

    /**
     * Returns the trypsin example as amino acid pattern.
     *
     * @return the trypsin example as amino acid pattern
     */
    public static AminoAcidPattern getTrypsinExample() {
        AminoAcidPattern example = new AminoAcidPattern();
        example.setTarget(0);
        ArrayList<Character> target = new ArrayList<>(2);
        target.add(AminoAcid.K.getSingleLetterCodeAsChar());
        target.add(AminoAcid.R.getSingleLetterCodeAsChar());
        example.setTargeted(0, target);
        ArrayList<Character> exclusion = new ArrayList<>(1);
        exclusion.add(AminoAcid.P.getSingleLetterCodeAsChar());
        example.setExcluded(1, exclusion);
        return example;
    }

    /**
     * Simple merger for two patterns.
     *
     * Example: this: target{0&gt;S} otherPattern: target{0&gt;T} result (this):
     * target{0&gt;S|T}
     *
     * @param otherPattern another pattern to be merged with this
     */
    public void merge(AminoAcidPattern otherPattern) {

        HashMap<Integer, ArrayList<Character>> otherInclusionMap = otherPattern.getAaTargeted();

        if (otherInclusionMap != null) {
            for (int i : otherInclusionMap.keySet()) {
                ArrayList<Character> otherAAs = otherPattern.getTargetedAA(i);
                if (residueTargeted == null) {
                    residueTargeted = new HashMap<>(otherInclusionMap.size());
                }
                ArrayList<Character> targetedAA = residueTargeted.get(i);
                if (targetedAA == null) {
                    residueTargeted.put(i, (ArrayList<Character>) otherAAs.clone());
                } else if (!otherAAs.isEmpty()) {
                    for (Character aa : otherAAs) {
                        if (!targetedAA.contains(aa)) {
                            targetedAA.add(aa);
                        }
                    }
                } else {
                    targetedAA.clear();
                }
                if (i + 1 > length) {
                    length = i + 1;
                }
            }
        }
    }

    /**
     * Appends another pattern at the end of this pattern.
     *
     * @param otherPattern the other pattern to append.
     */
    public void append(AminoAcidPattern otherPattern) {

        int patternLength = length();
        HashMap<Integer, ArrayList<Character>> otherTargetedMap = otherPattern.getAaTargeted();
        if (otherTargetedMap != null) {
            if (residueTargeted == null) {
                residueTargeted = new HashMap<>(otherTargetedMap.size());
            }
            for (int i : otherTargetedMap.keySet()) {
                int index = patternLength + i;
                residueTargeted.put(index, (ArrayList<Character>) otherTargetedMap.get(i).clone());
            }
        }

        length = patternLength + otherPattern.length();
    }

    /**
     * Convenience method merging two different patterns (see public void
     * merge(AminoAcidPattern otherPattern) for detailed information of the
     * merging procedure).
     *
     * @param pattern1 the first pattern
     * @param pattern2 the second pattern
     * @return a merged version of the two patterns
     */
    public static AminoAcidPattern merge(AminoAcidPattern pattern1, AminoAcidPattern pattern2) {
        AminoAcidPattern result = new AminoAcidPattern(pattern1);
        result.merge(pattern2);
        return result;
    }

    @Override
    public String toString() {
        return asStringBuilder().toString();
    }

    /**
     * Returns the sequence represented by this amino acid pattern in a new
     * string builder.
     *
     * @return the sequence represented by this amino acid pattern in a new
     * string builder
     */
    public StringBuilder asStringBuilder() {
        StringBuilder result = new StringBuilder(length());
        for (int i = 0; i < length(); i++) {
            if (getNTargetedAA(i) == 1) {
                result.append(getTargetedAA(i).get(0));
            } else {
                int nTargetedAas = getNTargetedAA(i);
                switch (nTargetedAas) {
                    case 0:
                        result.append("X");
                        break;
                    case 1:
                        result.append(getTargetedAA(i).get(0));
                        break;
                    default:
                        result.append("[");
                        for (Character aa : getTargetedAA(i)) {
                            result.append(aa);
                        }
                        result.append("]");
                        break;
                }
            }
        }
        return result;
    }

    /**
     * Returns the component of the amino acid pattern at the given index. 0 is
     * the first amino acid.
     *
     * @param index the index in the pattern. 0 is the first amino acid
     *
     * @return the component of the amino acid pattern at the given index
     */
    public String asSequence(int index) {
        return asStringBuilder().substring(index, index + 1);
    }

    /**
     * Adds a list of modifications to one of the amino acid pattern.
     *
     * @param localization the index of the amino acid residue site
     * @param ModificationSite valid amino acids for this site
     */
    public void addModificationSite(int localization, ArrayList<Character> ModificationSite) {
        if (residueTargeted == null) {
            residueTargeted = new HashMap<>(1);
        }
        residueTargeted.put(localization, ModificationSite);
    }

    /**
     * Returns all possible sequences which can be obtained from the targeted
     * amino acids. Missing amino acids will be denoted as 'X'. This does not
     * implement excluded amino acids.
     *
     * @return all possible sequences which can be obtained from the targeted
     * amino acids
     */
    public ArrayList<String> getAllPossibleSequences() {

        ArrayList<StringBuilder> stringBuilders = new ArrayList<>(1);
        int tempLength = length();

        for (int i = 0; i < length(); i++) {
            if (stringBuilders.isEmpty()) {
                if (residueTargeted != null) {
                    ArrayList<Character> aminoAcids = residueTargeted.get(i);
                    if (aminoAcids != null && !aminoAcids.isEmpty()) {
                        for (Character aminoAcid : aminoAcids) {
                            StringBuilder newBuilder = new StringBuilder(tempLength);
                            newBuilder.append(aminoAcid);
                            stringBuilders.add(newBuilder);
                        }
                    } else {
                        StringBuilder newBuilder = new StringBuilder(tempLength);
                        newBuilder.append('X');
                        stringBuilders.add(newBuilder);
                    }
                } else {
                    StringBuilder newBuilder = new StringBuilder(tempLength);
                    newBuilder.append('X');
                    stringBuilders.add(newBuilder);
                }
            } else {
                ArrayList<StringBuilder> newBuilders = new ArrayList<>(1);
                for (StringBuilder stringBuilder : stringBuilders) {
                    if (residueTargeted != null) {
                        ArrayList<Character> aminoAcids = residueTargeted.get(i);
                        if (aminoAcids != null && !aminoAcids.isEmpty()) {
                            for (Character aminoAcid : residueTargeted.get(i)) {
                                StringBuilder newBuilder = new StringBuilder(stringBuilder);
                                newBuilder.append(aminoAcid);
                                newBuilders.add(newBuilder);
                            }
                        } else {
                            StringBuilder newBuilder = new StringBuilder(stringBuilder);
                            newBuilder.append('X');
                            newBuilders.add(newBuilder);
                        }
                    } else {
                        StringBuilder newBuilder = new StringBuilder(stringBuilder);
                        newBuilder.append('X');
                        newBuilders.add(newBuilder);
                    }
                }
                stringBuilders = newBuilders;
            }
        }

        ArrayList<String> results = new ArrayList<>(stringBuilders.size());
        for (StringBuilder stringBuilder : stringBuilders) {
            results.add(stringBuilder.toString());
        }

        return results;
    }

    /**
     * Returns a sub pattern of the pattern.
     *
     * @param startIndex the start index, inclusive (0 is the first amino acid)
     * @param endIndex the end index, inclusive
     * @param updateTarget boolean indicating whether the target of the pattern
     * shall be updated. If yes it will be shifted by startIndex, simply copied
     * otherwise.
     *
     * @return a sub pattern
     */
    public AminoAcidPattern getSubPattern(int startIndex, int endIndex, boolean updateTarget) {

        AminoAcidPattern aminoAcidPattern = new AminoAcidPattern();

        if (residueTargeted != null) {
            for (int i : residueTargeted.keySet()) {
                if (i >= startIndex && i <= endIndex) {
                    ArrayList<Character> aminoAcids = (ArrayList<Character>) residueTargeted.get(i).clone();
                    aminoAcidPattern.setTargeted(i - startIndex, aminoAcids);
                }
            }
        }

        if (updateTarget) {
            aminoAcidPattern.setTarget(getTarget() - startIndex);
        } else {
            aminoAcidPattern.setTarget(getTarget());
        }

        return aminoAcidPattern;
    }

    /**
     * Returns a sub pattern of the pattern.
     *
     * @param startIndex the start index, inclusive (0 is the first amino acid)
     * @param updateTarget boolean indicating whether the target of the pattern
     * shall be updated. If yes it will be shifted by startIndex, simply copied
     * otherwise.
     *
     * @return a sub pattern
     */
    public AminoAcidPattern getSubPattern(int startIndex, boolean updateTarget) {
        return getSubPattern(startIndex, length(), updateTarget);
    }

    /**
     * Returns an amino acid pattern which is a reversed version of the current
     * pattern.
     *
     * @return an amino acid pattern which is a reversed version of the current
     * pattern
     */
    public AminoAcidPattern reverse() {

        AminoAcidPattern newPattern = new AminoAcidPattern();

        if (residueTargeted != null) {
            for (int i : residueTargeted.keySet()) {
                int reversed = length() - i - 1;
                newPattern.setTargeted(reversed, (ArrayList<Character>) residueTargeted.get(i).clone());
            }
        }

        /*
        if (target > -1) {
            newPattern.setTarget(length() - target - 1);
        }
         */
        return newPattern;
    }
}
