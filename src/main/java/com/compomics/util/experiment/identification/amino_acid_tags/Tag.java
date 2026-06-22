package com.compomics.util.experiment.identification.amino_acid_tags;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidPattern;
import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidSequence;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.utils.ModificationUtils;
import com.compomics.util.experiment.mass_spectrometry.utils.StandardMasses;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.Collectors;

/**
 * This class represents a sequence mass tag.
 *
 * @author Marc Vaudel
 */
public class Tag extends ExperimentObject {

    /**
     * Serial number for backward compatibility.
     */
    static final long serialVersionUID = 1625541843008045218L;
    /**
     * The content of the tag.
     */
    private ArrayList<TagComponent> content = new ArrayList<>(1);

    /**
     * Constructor for an empty tag.
     */
    public Tag() {

    }

    /**
     * Creates a new tag instance based on the given one.
     *
     * @param tag the reference tag
     */
    public Tag(Tag tag) {

        for (TagComponent tagComponent : tag.getContent()) {

            if (tagComponent instanceof MassGap) {

                MassGap massGap = (MassGap) tagComponent;
                addMassGap(massGap.getMass());

            } else if (tagComponent instanceof AminoAcidSequence) {

                AminoAcidSequence aminoAcidSequence = new AminoAcidSequence((AminoAcidSequence) tagComponent);
                addAminoAcidSequence(aminoAcidSequence);

            } else {

                throw new UnsupportedOperationException("Tag constructor not implemeted for tag component " + tagComponent.getClass() + ".");

            }
        }
    }

    /**
     * Constructor for a tag consisting of a sequence tag between two mass tags.
     *
     * @param nTermGap the N-term mass gap
     * @param sequenceTag the sequence tag with modifications
     * @param cTermGap the C-term mass gap
     */
    public Tag(double nTermGap, AminoAcidSequence sequenceTag, double cTermGap) {

        addMassGap(nTermGap);
        addAminoAcidSequence(sequenceTag);
        addMassGap(cTermGap);

    }

    /**
     * Returns the content of this tag as a list.
     *
     * @return the content of this tag as a list
     */
    public ArrayList<TagComponent> getContent() {

        return content;

    }

    /**
     * Sets the content for the given tag.
     *
     * @param content an array of tag components
     */
    public void setContent(ArrayList<TagComponent> content) {

        this.content = content;

    }

    /**
     * Adds a mass gap to the tag.
     *
     * @param massGap the value of the mass gap
     */
    public void addMassGap(double massGap) {

        if (massGap != 0) {

            content.add(new MassGap(massGap));

        }
    }

    /**
     * Adds a sequence of amino acids to the tag.
     *
     * @param aminoAcidSequence the amino acid sequence with modifications
     */
    public void addAminoAcidSequence(AminoAcidSequence aminoAcidSequence) {

        if (aminoAcidSequence.length() > 0) {

            if (!content.isEmpty()) {

                TagComponent lastComponent = content.get(content.size() - 1);

                if (lastComponent instanceof AminoAcidSequence) {

                    AminoAcidSequence sequence = (AminoAcidSequence) lastComponent;
                    sequence.appendCTerm(aminoAcidSequence);
                    return;

                }
            }

            content.add(aminoAcidSequence);

        }
    }

    /**
     * Returns the tag as intelligible sequence of tag components. For example
     * amino acid tags and mass gaps: &lt;115.2&gt;TAG&lt;110.5&gt;.
     *
     * @return The tag as intelligible sequence for display.
     */
    public String asSequence() {

        return content.stream()
                .map(TagComponent::asSequence)
                .collect(Collectors.joining());

    }

    /**
     * Returns the longest amino acid sequence contained in this tag.
     *
     * @return the longest amino acid sequence contained in this tag
     */
    public String getLongestAminoAcidSequence() {

        String result = "";
        AminoAcidSequence lastAminoAcidSequence = null;

        for (TagComponent tagComponent : content) {

            if (tagComponent instanceof MassGap) {

                if (lastAminoAcidSequence != null && lastAminoAcidSequence.length() > result.length()) {

                    result = lastAminoAcidSequence.asSequence();

                }

                lastAminoAcidSequence = null;

            } else if (tagComponent instanceof AminoAcidSequence) {

                AminoAcidSequence currentSequence = (AminoAcidSequence) tagComponent;

                if (lastAminoAcidSequence != null) {

                    lastAminoAcidSequence.appendCTerm(currentSequence);

                }

            } else {

                throw new UnsupportedOperationException(
                        "Longest amino acid sequence not implemented for tag component "
                        + tagComponent.getClass()
                        + "."
                );

            }
        }

        if (lastAminoAcidSequence != null && lastAminoAcidSequence.length() > result.length()) {

            result = lastAminoAcidSequence.asSequence();

        }

        return result;

    }

    /**
     * Returns the mass of the tag.
     *
     * @return the mass of the tag
     */
    public double getMass() {

        return StandardMasses.h2o.mass
                + content.stream()
                        .mapToDouble(tagComponent -> tagComponent.getMass())
                        .sum();

    }

    /**
     * Returns the theoretic mass of the tag, eventually without terminal gaps.
     *
     * @param includeCTermGap if true the C-terminal gap will be added if
     * present
     * @param includeNTermGap if true the N-terminal gap will be added if
     * present
     *
     * @return the theoretic mass of the tag
     */
    public double getMass(boolean includeCTermGap, boolean includeNTermGap) {

        double mass = getMass();

        if (!includeCTermGap) {

            mass -= getCTerminalGap();

        }

        if (!includeNTermGap) {

            mass -= getNTerminalGap();

        }

        return mass;

    }

    /**
     * Returns the N-terminal gap of the tag.
     *
     * @return the N-terminal gap of the tag
     */
    public double getNTerminalGap() {

        if (content.isEmpty()) {

            return 0.0;

        }

        TagComponent tagComponent = content.get(0);

        if (tagComponent instanceof MassGap) {

            return tagComponent.getMass();

        } else {

            return 0.0;

        }
    }

    /**
     * Returns the C-terminal gap of the tag.
     *
     * @return the C-terminal gap of the tag
     */
    public double getCTerminalGap() {

        if (content.isEmpty()) {

            return 0.0;

        }

        TagComponent tagComponent = content.get(content.size() - 1);

        if (tagComponent instanceof MassGap) {

            return tagComponent.getMass();

        } else {

            return 0.0;

        }
    }

    /**
     * Returns the modified sequence as an tagged string with potential
     * modification sites color coded or with modification tags, e.g,
     * &lt;mox&gt;. /!\ this method will work only if the modifications found in
     * the tag components are in the factory. /!\ This method uses the
     * modifications as set in the modification matches of this peptide and
     * displays all of them.
     *
     * @param modificationProfile the modification profile of the search
     * @param useHtmlColorCoding if true, color coded HTML is used, otherwise
     * modification tags, e.g, &lt;mox&gt;, are used
     * @param includeHtmlStartEndTags if true, start and end HTML tags are added
     * @param useShortName if true, the short names are used in the tags
     * @param includeTerminalGaps if true, the terminal gaps will be displayed
     * on the sequence
     * @param modificationsSequenceMatchingParameters the sequence matching
     * parameters to use for modifications
     * @param displayedModifications the modifications to display
     *
     * @return the modified sequence as a tagged string
     */
    public String getTaggedModifiedSequence(
            ModificationParameters modificationProfile,
            boolean useHtmlColorCoding,
            boolean includeHtmlStartEndTags,
            boolean useShortName,
            boolean includeTerminalGaps,
            SequenceMatchingParameters modificationsSequenceMatchingParameters,
            HashSet<String> displayedModifications
    ) {

        return getTaggedModifiedSequence(
                modificationProfile,
                this,
                useHtmlColorCoding,
                includeHtmlStartEndTags,
                useShortName,
                includeTerminalGaps,
                modificationsSequenceMatchingParameters,
                displayedModifications
        );

    }

    public String getTaggedModifiedSequence(
            ModificationParameters modificationProfile,
            boolean useHtmlColorCoding,
            boolean includeHtmlStartEndTags,
            boolean useShortName,
            boolean includeTerminalGaps
    ) {
        return getTaggedModifiedSequence(
                modificationProfile,
                useHtmlColorCoding,
                includeHtmlStartEndTags,
                useShortName,
                includeTerminalGaps,
                SequenceMatchingParameters.DEFAULT_STRING_MATCHING,
                null
        );
    }

    public String getTaggedModifiedSequence(
            ModificationParameters modificationProfile,
            boolean useHtmlColorCoding,
            boolean includeHtmlStartEndTags,
            boolean useShortName,
            boolean includeTerminalGaps,
            boolean excludeAllFixedPtms
    ) {
        return getTaggedModifiedSequence(
                modificationProfile,
                useHtmlColorCoding,
                includeHtmlStartEndTags,
                useShortName,
                includeTerminalGaps,
                SequenceMatchingParameters.DEFAULT_STRING_MATCHING,
                null
        );
    }

    /**
     * Returns the modified sequence as an tagged string with potential
     * modification sites color coded or with modification tags, e.g,
     * &lt;mox&gt;. /!\ This method will work only if the modification found in
     * the peptide are in the factory. /!\ This method uses the modifications as
     * set in the modification matches of this peptide and displays all of them.
     *
     * @param modificationParameters the modification profile of the search
     * @param tag the tag
     * @param includeHtmlStartEndTags if true, start and end HTML tags are added
     * @param useHtmlColorCoding if true, color coded HTML is used, otherwise
     * modification tags, e.g, &lt;mox&gt;, are used
     * @param useShortName if true, the short names are used in the tags
     * @return the tagged modified sequence as a string
     * @param includeTerminalGaps if true, the terminal gaps will be displayed
     * on the sequence
     * @param modificationsSequenceMatchingParameters the sequence matching
     * parameters to use for modifications
     * @param displayedModifications the modifications to display
     */
    public static String getTaggedModifiedSequence(
            ModificationParameters modificationParameters,
            Tag tag,
            boolean useHtmlColorCoding,
            boolean includeHtmlStartEndTags,
            boolean useShortName,
            boolean includeTerminalGaps,
            SequenceMatchingParameters modificationsSequenceMatchingParameters,
            HashSet<String> displayedModifications
    ) {

        StringBuilder modifiedSequence = new StringBuilder();

        if (useHtmlColorCoding && includeHtmlStartEndTags) {
            modifiedSequence.append("<html>");
        }

        modifiedSequence.append(tag.getNTerminal(includeTerminalGaps, modificationParameters, modificationsSequenceMatchingParameters)).append('-');

        for (int i = 0; i < tag.getContent().size(); i++) {

            TagComponent tagComponent = tag.getContent().get(i);

            if (tagComponent instanceof AminoAcidSequence) {

                AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) tagComponent;

                String[] variableModifications = aminoAcidSequence.getIndexedVariableModifications();

                // remove the hidden modifications
                for (int j = 0; j < variableModifications.length; j++) { // @TODO: possible to do this with streams?
                    String tempMod = variableModifications[j];
                    if (tempMod != null && displayedModifications != null && !displayedModifications.contains(tempMod)) {
                        variableModifications[j] = null;
                    }
                }

                String[] fixedModifications = aminoAcidSequence.getFixedModifications(i == 0, i == tag.getContent().size() - 1, modificationParameters, modificationsSequenceMatchingParameters);

                // remove the hidden modifications
                for (int j = 0; j < fixedModifications.length; j++) { // @TODO: possible to do this with streams?
                    String tempMod = fixedModifications[j];
                    if (tempMod != null && displayedModifications != null && !displayedModifications.contains(tempMod)) {
                        fixedModifications[j] = null;
                    }
                }

                modifiedSequence.append(ModificationUtils.getTaggedModifiedSequence(modificationParameters, aminoAcidSequence.getSequence(), variableModifications, null, null, fixedModifications, useHtmlColorCoding, useShortName));

            } else if (tagComponent instanceof MassGap) {

                if (includeTerminalGaps || i > 0 && i < tag.getContent().size() - 1) {
                    modifiedSequence.append(tagComponent.asSequence());
                }

            } else {
                throw new UnsupportedOperationException("Tagged sequence not implemented for tag component " + tagComponent.getClass() + ".");
            }
        }

        modifiedSequence.append('-').append(tag.getCTerminal(includeTerminalGaps, modificationParameters, modificationsSequenceMatchingParameters));

        if (useHtmlColorCoding && includeHtmlStartEndTags) {
            modifiedSequence.append("</html>");
        }

        return modifiedSequence.toString();
    }

    /**
     * Returns the N-terminal tag of this tag as a string for sequence display.
     *
     * @param includeTerminalGaps indicates whether mass gaps shall be included
     * @param modificationParameters the modification profile of the search
     * @param modificationsSequenceMatchingParameters the sequence matching
     * parameters to use for modifications
     *
     * @return the N-terminal tag of this tag
     */
    public String getNTerminal(
            boolean includeTerminalGaps,
            ModificationParameters modificationParameters,
            SequenceMatchingParameters modificationsSequenceMatchingParameters
    ) {

        if (content.isEmpty()) {
            return "";
        }

        TagComponent firstComponent = content.get(0);

        if (firstComponent instanceof MassGap) {

            return includeTerminalGaps ? firstComponent.asSequence() : Double.toString(Util.roundDouble(firstComponent.getMass(), 2));

        } else if (firstComponent instanceof AminoAcidSequence) {

            AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) firstComponent;

            String nTermMod = aminoAcidSequence.getFixedModifications(true, false, modificationParameters, modificationsSequenceMatchingParameters)[0];

            if (nTermMod != null) {
                return Util.replaceAll(nTermMod, '-', ' ');
            }

            nTermMod = aminoAcidSequence.getIndexedVariableModifications()[0];

            if (nTermMod != null) {
                return Util.replaceAll(nTermMod, '-', ' ');
            }

            return "NH2";

        } else {
            throw new UnsupportedOperationException("N-terminal tag not implemented for tag component " + firstComponent.getClass() + ".");
        }
    }

    /**
     * Returns the C-terminal tag of this tag as a string for sequence display.
     *
     * @param includeTerminalGaps indicates whether mass gaps shall be included
     * @param modificationParameters the modification profile of the search
     * @param modificationsSequenceMatchingParameters the sequence matching
     * parameters to use for modifications
     *
     * @return the C-terminal tag of this tag
     */
    public String getCTerminal(
            boolean includeTerminalGaps,
            ModificationParameters modificationParameters,
            SequenceMatchingParameters modificationsSequenceMatchingParameters
    ) {

        if (content.isEmpty()) {

            return "";

        }

        TagComponent lastComponent = content.get(content.size() - 1);

        if (lastComponent instanceof MassGap) {

            return includeTerminalGaps ? lastComponent.asSequence() : Double.toString(Util.roundDouble(lastComponent.getMass(), 2));

        } else if (lastComponent instanceof AminoAcidSequence) {

            AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) lastComponent;

            String cTermMod = aminoAcidSequence.getFixedModifications(true, false, modificationParameters, modificationsSequenceMatchingParameters)[aminoAcidSequence.length() + 1];

            if (cTermMod != null) {
                return Util.replaceAll(cTermMod, '-', ' ');
            }

            cTermMod = aminoAcidSequence.getIndexedVariableModifications()[aminoAcidSequence.length() + 1];

            if (cTermMod != null) {
                return Util.replaceAll(cTermMod, '-', ' ');
            }

            return "COOH";

        } else {
            throw new UnsupportedOperationException("C-terminal tag not implemented for tag component " + lastComponent.getClass() + ".");
        }
    }

    /**
     * Returns the amino acid length of the tag when mass gaps are considered
     * like one amino acid
     *
     * @return the amino acid length of the tag
     */
    public int getLengthInAminoAcid() {

        int length = 0;

        for (TagComponent tagComponent : content) {

            if (tagComponent instanceof AminoAcidSequence) {

                AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) tagComponent;
                length += aminoAcidSequence.length();

            } else if (tagComponent instanceof MassGap) {

                length++;

            } else {

                throw new UnsupportedOperationException("Tag length in amino acid not implemented for tag component " + tagComponent.getClass() + ".");

            }
        }

        return length;

    }

    /**
     * Returns the potential modification sites as an ordered list of string. 1
     * is the first amino acid. An empty list is returned if no possibility was
     * found. This method does not account for protein terminal modifications.
     *
     * @param modification the modification considered
     * @param modificationSequenceMatchingPreferences the sequence matching
     * preferences for the modification to amino acid sequence mapping
     *
     * @return a list of potential modification sites
     */
    public ArrayList<Integer> getPotentialModificationSites(
            Modification modification,
            SequenceMatchingParameters modificationSequenceMatchingPreferences
    ) {

        ArrayList<Integer> possibleSites = new ArrayList<>();
        AminoAcidPattern modificationPattern = modification.getPattern();
        int patternLength = modificationPattern.length();

        switch (modification.getModificationType()) {

            case modaa:

                int offset = 0;

                for (TagComponent tagComponent : content) {

                    if (tagComponent instanceof AminoAcidSequence) {

                        AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) tagComponent;

                        for (int i : modificationPattern.getIndexes(aminoAcidSequence.getSequence(), modificationSequenceMatchingPreferences)) {

                            possibleSites.add(i + offset);

                        }

                        offset += aminoAcidSequence.length();

                    } else {

                        offset++;

                    }
                }

                return possibleSites;

            case modc_protein:
            case modc_peptide:

                possibleSites.add(patternLength);

                return possibleSites;

            case modn_peptide:
            case modn_protein:

                possibleSites.add(1);
                return possibleSites;

            case modcaa_peptide:
            case modcaa_protein:

                if (content.isEmpty()) {

                    return new ArrayList<>(0);

                }

                TagComponent component = content.get(content.size() - 1);

                if (component instanceof AminoAcidSequence) {

                    AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) component;

                    if (modificationPattern.isEnding(aminoAcidSequence.getSequence(), modificationSequenceMatchingPreferences)) {

                        possibleSites.add(patternLength);

                    }

                } else if (component instanceof MassGap) {

                    possibleSites.add(patternLength);

                } else {

                    throw new UnsupportedOperationException("Possible modifications not implemnted for tag component " + component.getClass() + ".");

                }

                return possibleSites;

            case modnaa_peptide:
            case modnaa_protein:

                if (content.isEmpty()) {

                    return new ArrayList<>(0);

                }

                component = content.get(0);

                if (component instanceof AminoAcidSequence) {

                    AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) component;

                    if (modificationPattern.isStarting(aminoAcidSequence.getSequence(), modificationSequenceMatchingPreferences)) {

                        possibleSites.add(patternLength);

                    }

                } else if (component instanceof MassGap) {

                    possibleSites.add(patternLength);

                } else {

                    throw new UnsupportedOperationException("Possible modifications not implemnted for tag component " + component.getClass() + ".");

                }

                return possibleSites;

            default:

                throw new IllegalArgumentException("Modification type " + modification.getModificationType() + " not recognized.");

        }
    }

    /**
     * Indicates whether this tag is the same as another tag. Note: this method
     * accounts for modification localization.
     *
     * @param anotherTag another tag
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return a boolean indicating whether the tag is the same as another
     */
    public boolean isSameAs(
            Tag anotherTag, 
            SequenceMatchingParameters sequenceMatchingPreferences
    ) {

        if (content.size() != anotherTag.getContent().size()) {

            return false;

        }

        for (int i = 0; i < content.size(); i++) {

            TagComponent component1 = content.get(i);
            TagComponent component2 = anotherTag.getContent().get(i);

            if (!component1.isSameAs(component2, sequenceMatchingPreferences)) {

                return false;

            }
        }

        return true;

    }

    /**
     * Indicates whether this tag is the same as another tag without accounting
     * for modification localization.
     *
     * @param anotherTag another tag
     * @param sequenceMatchingPreferences the sequence matching preferences
     *
     * @return a boolean indicating whether the tag is the same as another
     */
    public boolean isSameSequenceAndModificationStatusAs(
            Tag anotherTag, 
            SequenceMatchingParameters sequenceMatchingPreferences
    ) {

        if (content.size() != anotherTag.getContent().size()) {

            return false;

        }

        for (int i = 0; i < content.size(); i++) {

            TagComponent component1 = content.get(i);
            TagComponent component2 = anotherTag.getContent().get(i);

            if (!component1.isSameSequenceAndModificationStatusAs(component2, sequenceMatchingPreferences)) {

                return false;

            }
        }

        return true;

    }

    /**
     * Returns the tag modifications as a string.
     *
     * @param tag the tag
     *
     * @return the peptide modifications as a string
     */
    public static String getTagModificationsAsString(Tag tag) {

        HashMap<String, ArrayList<Integer>> modMap = new HashMap<>();
        int offset = 0;

        for (TagComponent tagComponent : tag.getContent()) {

            if (tagComponent instanceof MassGap) {

                offset++;

            } else if (tagComponent instanceof AminoAcidSequence) {

                AminoAcidSequence aminoAcidSequence = (AminoAcidSequence) tagComponent;

                ModificationMatch[] variableModifications = aminoAcidSequence.getVariableModifications();

                for (ModificationMatch modificationMatch : variableModifications) {

                    if (!modMap.containsKey(modificationMatch.getModification())) {

                        modMap.put(modificationMatch.getModification(), new ArrayList<>(1));

                    }

                    modMap.get(modificationMatch.getModification()).add(modificationMatch.getSite() + offset);

                }

                offset += aminoAcidSequence.length();

            } else {

                throw new UnsupportedOperationException("Modification summary not implemented for TagComponent " + tagComponent.getClass() + ".");

            }
        }

        return modMap.entrySet().stream()
                .map(entry -> getModificationLine(entry.getKey(), entry.getValue()))
                .collect(Collectors.joining(","));
    }

    /**
     * Returns a description line for a modification and its sites.
     *
     * @param modName the name of the modification
     * @param sites the list of sites
     *
     * @return a description line for a modification and its sites
     */
    private static String getModificationLine(String modName, ArrayList<Integer> sites) {

        String sitesString = sites.stream()
                .map(site -> site.toString())
                .collect(Collectors.joining(","));

        StringBuilder modLine = new StringBuilder(modName.length() + sitesString.length() + 2);
        modLine.append(modName).append("(").append(sitesString).append(")");

        return modLine.toString();
    }

    /**
     * Returns a new tag instance which is a reversed version of the current
     * tag.
     *
     * @param yIon indicates whether the tag is based on y ions
     *
     * @return a new tag instance which is a reversed version of the current tag
     */
    public Tag reverse(boolean yIon) {

        Tag newTag = new Tag();

        for (int i = content.size() - 1; i >= 0; i--) {

            TagComponent tagComponent = content.get(i);

            if (tagComponent instanceof MassGap) {

                double mass = tagComponent.getMass();

                if (i == content.size() - 1) {

                    if (yIon) {

                        mass += StandardMasses.h2o.mass;

                    } else {

                        mass -= StandardMasses.h2o.mass;

                    }

                } else if (i == 0) {

                    if (yIon) {

                        mass -= StandardMasses.h2o.mass;

                    } else {

                        mass += StandardMasses.h2o.mass;

                    }
                }

                newTag.addMassGap(mass);

            } else if (tagComponent instanceof AminoAcidSequence) {

                newTag.addAminoAcidSequence(((AminoAcidSequence) tagComponent).reverse());

            } else {

                throw new UnsupportedOperationException("Reverse method not implemented for tag component " + tagComponent.getClass() + ".");

            }
        }

        return newTag;

    }

    /**
     * Indicates whether the tag can be reversed (ie if termini are mass gaps
     * with mass &ge; water).
     *
     * @return whether the tag can be reversed
     */
    public boolean canReverse() {

        TagComponent terminalComponent = content.get(0);

        if (terminalComponent instanceof MassGap) {

            MassGap terminalGap = (MassGap) terminalComponent;

            if (terminalGap.getMass() >= StandardMasses.h2o.mass) {

                terminalComponent = content.get(content.size() - 1);

                if (terminalComponent instanceof MassGap) {

                    terminalGap = (MassGap) terminalComponent;

                    if (terminalGap.getMass() >= StandardMasses.h2o.mass) {

                        return true;

                    }
                }
            }
        }

        return false;

    }

    @Override
    public String toString() {
        return asSequence();
    }
}
