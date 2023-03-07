package com.compomics.util.experiment.io.identification;

import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.mass_spectrometry.SpectrumProvider;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.waiting.WaitingHandler;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import org.xmlpull.v1.XmlPullParserException;

/**
 * Generic interface for the parser of a file containing PSMs.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public interface IdfileReader extends AutoCloseable {

    /**
     * Returns the names and versions of the software used to generate the
     * identification file in a map, e.g., Mascot &gt; (2.2 and 2.3) and
     * X!Tandem &gt; Sledgehammer (2013.09.01.1). Null if not known.
     *
     * @return the version of the software used to generate the identification
     * file, null if not known
     */
    public HashMap<String, ArrayList<String>> getSoftwareVersions();

    /**
     * Returns the extension of the file for which this IdfileReader can be
     * used.
     *
     * @return String with the extension (taken to make up the end of the
     * filename) of the file that this IdfileReader can read.
     */
    public String getExtension();

    @Override
    public void close() throws IOException;

    /**
     * Retrieves all the spectrum matches from an identification file as a list
     * of spectrum matches, one spectrum match per spectrum. It is very
     * important to close the file reader after creation. Using this method
     * secondary maps are not filled.
     *
     * @param spectrumProvider A spectrum provider with the spectra of the file loaded.
     * @param waitingHandler The waiting handler displaying the progress (can be
     * null). The secondary progress methods will be called.
     * @param searchParameters The search parameters.
     *
     * @return a list of spectrum matches
     *
     * @throws IOException if an IOException occurs
     * @throws SQLException if an SQLException occurs
     * @throws ClassNotFoundException if an\ ClassNotFoundException occurs
     * @throws InterruptedException if an InterruptedException occurs
     * @throws JAXBException if a JAXBException occurs
     * @throws XmlPullParserException if an XmlPullParserException occurs
     * @throws XMLStreamException if an XMLStreamException occurs
     */
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler,
            SearchParameters searchParameters
    ) 
            throws IOException, SQLException, ClassNotFoundException, InterruptedException, JAXBException, XmlPullParserException, XMLStreamException;

    /**
     * Retrieves all the spectrum matches from an identification file as a list
     * of spectrum matches, one spectrum match per spectrum.It is very important
     * to close the file reader after creation. Secondary peptide and tag maps
     * are filled according to the file content and the sequence matching
     * preferences. If the sequence matching preferences are null, the maps are
     * not filled.
     *
     * @param spectrumProvider A spectrum provider with the spectra of the file loaded.
     * @param waitingHandler The waiting handler displaying the progress (can be
     * null). The secondary progress methods will be called.
     * @param searchParameters The search parameters.
     * @param sequenceMatchingPreferences The sequence matching preferences to
     * use for the creation of the secondary maps.
     * @param expandAaCombinations If true, a peptide assumption (not
     * implemented for tag assumptions) will be created for all possible amino
     * acid combination for peptide sequences containing an ambiguity like an X.
     *
     * @return the spectrum matches
     *
     * @throws IOException if an IOException occurs
     * @throws SQLException if an SQLException occurs
     * @throws ClassNotFoundException if an\ ClassNotFoundException occurs
     * @throws InterruptedException if an InterruptedException occurs
     * @throws JAXBException if a JAXBException occurs
     * @throws XmlPullParserException if an XmlPullParserException occurs
     * @throws XMLStreamException if an XMLStreamException occurs
     */
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler,
            SearchParameters searchParameters,
            SequenceMatchingParameters sequenceMatchingPreferences,
            boolean expandAaCombinations
    ) 
            throws IOException, SQLException, ClassNotFoundException, InterruptedException, JAXBException, XmlPullParserException, XmlPullParserException, XMLStreamException;

    /**
     * Returns a boolean indicating whether the file contains de novo results as
     * tags.
     *
     * @return a boolean indicating whether the file contains de novo results as
     * tags
     */
    public boolean hasDeNovoTags();
}
