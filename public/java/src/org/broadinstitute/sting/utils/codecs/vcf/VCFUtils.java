/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.vcf;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * A set of static utility methods for common operations on VCF files/records.
 */
public class VCFUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private VCFUtils() { }

    public static Map<String, VCFHeader> getVCFHeadersFromRods(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        Map<String, VCFHeader> data = new HashMap<String, VCFHeader>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if it's not in our list
            if ( rodNames != null && !rodNames.contains(source.getName()) )
                continue;

            if ( source.getHeader() != null && source.getHeader() instanceof VCFHeader )
                data.put(source.getName(), (VCFHeader)source.getHeader());
        }

        return data;
    }

    public static Map<String,VCFHeader> getVCFHeadersFromRodPrefix(GenomeAnalysisEngine toolkit,String prefix) {
        Map<String, VCFHeader> data = new HashMap<String, VCFHeader>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if lacks the prefix
            if ( ! source.getName().startsWith(prefix) )
                continue;

            if ( source.getHeader() != null && source.getHeader() instanceof VCFHeader )
                data.put(source.getName(), (VCFHeader)source.getHeader());
        }

        return data;
    }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit) {
        return getHeaderFields(toolkit, null);
    }

    /**
     * Gets the header fields from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     * @param rodNames   names of rods to use, or null if we should use all possible ones
     *
     * @return a set of all fields
     */
    public static Set<VCFHeaderLine> getHeaderFields(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {

        // keep a map of sample name to occurrences encountered
        TreeSet<VCFHeaderLine> fields = new TreeSet<VCFHeaderLine>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            // ignore the rod if it's not in our list
            if ( rodNames != null && !rodNames.contains(source.getName()) )
                continue;

            if ( source.getRecordType().equals(VariantContext.class)) {
                VCFHeader header = (VCFHeader)source.getHeader();
                if ( header != null )
                    fields.addAll(header.getMetaData());
            }
        }

        return fields;
    }

    /** Only displays a warning if a logger is provided and an identical warning hasn't been already issued */
    private static final class HeaderConflictWarner {
        Logger logger;
        Set<String> alreadyIssued = new HashSet<String>();

        private HeaderConflictWarner(final Logger logger) {
            this.logger = logger;
        }

        public void warn(final VCFHeaderLine line, final String msg) {
            if ( logger != null && ! alreadyIssued.contains(line.getKey()) ) {
                alreadyIssued.add(line.getKey());
                logger.warn(msg);
            }
        }
    }

    public static Set<VCFHeaderLine> smartMergeHeaders(Collection<VCFHeader> headers, Logger logger) throws IllegalStateException {
        HashMap<String, VCFHeaderLine> map = new HashMap<String, VCFHeaderLine>(); // from KEY.NAME -> line
        HeaderConflictWarner conflictWarner = new HeaderConflictWarner(logger);

        // todo -- needs to remove all version headers from sources and add its own VCF version line
        for ( VCFHeader source : headers ) {
            //System.out.printf("Merging in header %s%n", source);
            for ( VCFHeaderLine line : source.getMetaData()) {
                String key = line.getKey();

                if ( line instanceof VCFNamedHeaderLine)
                    key = key + "" + ((VCFNamedHeaderLine) line).getName();

                if ( map.containsKey(key) ) {
                    VCFHeaderLine other = map.get(key);
                    if ( line.equals(other) )
                        continue;
                    else if ( ! line.getClass().equals(other.getClass()) )
                        throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    else if ( line instanceof VCFFilterHeaderLine) {
                        String lineName = ((VCFFilterHeaderLine) line).getName();
                        String otherName = ((VCFFilterHeaderLine) other).getName();
                        if ( ! lineName.equals(otherName) )
                            throw new IllegalStateException("Incompatible header types: " + line + " " + other );
                    } else if ( line instanceof VCFCompoundHeaderLine ) {
                        VCFCompoundHeaderLine compLine = (VCFCompoundHeaderLine)line;
                        VCFCompoundHeaderLine compOther = (VCFCompoundHeaderLine)other;

                        // if the names are the same, but the values are different, we need to quit
                        if (! (compLine).equalsExcludingDescription(compOther) ) {
                            if ( compLine.getType().equals(compOther.getType()) ) {
                                // The Number entry is an Integer that describes the number of values that can be
                                // included with the INFO field. For example, if the INFO field contains a single
                                // number, then this value should be 1. However, if the INFO field describes a pair
                                // of numbers, then this value should be 2 and so on. If the number of possible
                                // values varies, is unknown, or is unbounded, then this value should be '.'.
                                conflictWarner.warn(line, "Promoting header field Number to . due to number differences in header lines: " + line + " " + other);
                                compOther.setNumberToUnbounded();
                            } else if ( compLine.getType() == VCFHeaderLineType.Integer && compOther.getType() == VCFHeaderLineType.Float ) {
                                // promote key to Float
                                conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                                map.put(key, compOther);
                            } else if ( compLine.getType() == VCFHeaderLineType.Float && compOther.getType() == VCFHeaderLineType.Integer ) {
                                // promote key to Float
                                conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                            } else {
                                throw new IllegalStateException("Incompatible header types, collision between these two types: " + line + " " + other );
                            }
                        }
                        if ( ! compLine.getDescription().equals(compOther.getDescription()) )
                            conflictWarner.warn(line, "Allowing unequal description fields through: keeping " + compOther + " excluding " + compLine);
                    } else {
                        // we are not equal, but we're not anything special either
                        conflictWarner.warn(line, "Ignoring header line already in map: this header line = " + line + " already present header = " + other);
                    }
                } else {
                    map.put(key, line);
                    //System.out.printf("Adding header line %s%n", line);
                }
            }
        }

        return new HashSet<VCFHeaderLine>(map.values());
    }

    /**
     * Merges what is believed to be a CompoundHeaderLine(from INFO or FORMAT) into pre-existing entry
     *
     * @param map       the map to store the VCFHeaderLine in
     * @param key       the key formatted to be based on the id
     * @param line      the line that needs to be merged in
     * @param logger    logger to hold warning
     */
    private static void mergeCompoundHeaderLine(Map<String, VCFHeaderLine> map, String key, VCFHeaderLine line, Logger logger) throws IllegalStateException {
        HeaderConflictWarner conflictWarner = new HeaderConflictWarner(logger);
        VCFHeaderLine other = map.get(key);
        if(line.equals(other))
            return;
        if(!line.getClass().equals(other.getClass()))
            throw new IllegalStateException("Incompatible header types: " + line + " " + other);
        if(!(line instanceof VCFCompoundHeaderLine))
            throw new IllegalStateException("Expected VCFCompoundHeaderLine: " + line);

        VCFCompoundHeaderLine compLine = (VCFCompoundHeaderLine)line;
        VCFCompoundHeaderLine compOther = (VCFCompoundHeaderLine)other;

        if(compLine.getCount() != compOther.getCount()) {
            if ( logger != null ) conflictWarner.warn(line, "Promoting header field Number to . due to number differences in header lines: " + line + " " + other);
            compOther.setNumberToUnbounded();
        }

        if (compLine.getType() != (compOther.getType())) {
            if ( compLine.getType() == VCFHeaderLineType.Integer && compOther.getType() == VCFHeaderLineType.Float ) {
                if ( logger != null ) conflictWarner.warn(line, "Promoting Integer to Float in header: " + compLine);
            } else if ( compLine.getType() == VCFHeaderLineType.Float && compOther.getType() == VCFHeaderLineType.Integer ) {
                if ( logger != null ) conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
                compOther.promoteIntToFloat();
            } else {
                throw new IllegalStateException("Incompatible header types, collision between these two types: " + line + " " + other );
            }
        }

        if (!compLine.getDescription().equals(compOther.getDescription()))
            if ( logger != null ) conflictWarner.warn(line, "Allowing unequal description fields through: keeping " + compOther + " excluding " + compLine);

    }

    /**
     * Parses the ID out of the value of a "SAMPLE" VCFHeaderLine
     *
     * @param sampleString  the value of a "SAMPLE" VCFHeaderLine
     * @return              the ID field
     */

    private static String parseID(String sampleString) {
        int equalPos = sampleString.indexOf('=');
        int commaPos = sampleString.indexOf(',');
        return sampleString.substring(equalPos+1,commaPos);
    }

    /**
     * Inserts center with a . immediately after the first fields value (should be ID)
     *
     * @param sampleString  the value of a "SAMPLE" VCFHeaderLine
     * @param center        the name of the center
     * @return              a string such that ID=name becomes ID=name.center
     */
    private static String insertCenter(String sampleString, String center){
        return sampleString.replaceFirst(",", "." + center+ ",");
    }

    /**
     * Merges headers in a way more compliant with TCGA specifications, based on smartMergeHeaders
     *
     * This, and the associated functions above are probably too special cased and could use refactoring for
     * readability and reusability
     *
     * @param vcfRods   full rod list which includes the rod names
     * @param logger    logger to hold warnings
     * @return          set containing all the headers of the merged VCF
     * @throws IllegalStateException    generic exception if anything goes wrong
     */
    public static Set<VCFHeaderLine> tcgaMergeHeaders(Map<String, VCFHeader> vcfRods, Logger logger) throws IllegalStateException {
        HeaderConflictWarner conflictWarner = new HeaderConflictWarner(logger);
        HashMap<String, VCFHeaderLine> map = new HashMap<String, VCFHeaderLine>();
        HashMap<String, VCFHeaderLine> infoMap = new HashMap<String, VCFHeaderLine>();
        HashMap<String, VCFHeaderLine> filterMap = new HashMap<String, VCFHeaderLine>();
        HashMap<String, VCFHeaderLine> formatMap = new HashMap<String, VCFHeaderLine>();
        HashMap<String, VCFHeaderLine> sampleMap = new HashMap<String, VCFHeaderLine>();

        for(Map.Entry<String, VCFHeader> vcfRod : vcfRods.entrySet()) {
            String name = vcfRod.getKey();
            VCFHeader source = vcfRod.getValue();

            for (VCFHeaderLine line : source.getMetaData()) {
                String key = line.getKey();

                if(key.equals("center")) {
                    if(map.containsKey("center")) {
                        String newValue = map.get("center").getValue() + "," + line.getValue();
                        VCFHeaderLine newVCFHeaderLine = new VCFHeaderLine("center", newValue);
                        map.put("center", newVCFHeaderLine); // replace the old "center" with the new one
                    } else {
                        map.put("center", line);
                    }
                } else if(key.equals("INFO")) {
                    if(line instanceof VCFNamedHeaderLine) {
                        String infoKey = ((VCFNamedHeaderLine) line).getName();
                        if(infoMap.containsKey(infoKey))
                            mergeCompoundHeaderLine(infoMap, infoKey, line, logger);
                        else
                            infoMap.put(infoKey, line);
                    } else
                        throw new IllegalStateException("Incompatible header type, not VCFNamedHeaderLine: " + line);
                } else if(key.equals("FILTER")) {
                    if(line instanceof VCFFilterHeaderLine) {
                        //create a new VCFFilterHeaderLine with the annotated name
                        VCFFilterHeaderLine filterLine = (VCFFilterHeaderLine) line;
                        VCFFilterHeaderLine newLine = new VCFFilterHeaderLine(filterLine.getName() + "." + name, filterLine.getDescription());

                        String filterKey = newLine.getName();
                        if(filterMap.containsKey(filterKey)) {
                            VCFHeaderLine other = filterMap.get(filterKey);
                            String lineName = newLine.getName();
                            String otherName = ((VCFFilterHeaderLine) other).getName();
                            if ( !lineName.equals(otherName) )
                                throw new IllegalStateException("Incompatible header types: " + newLine + " " + other );
                        } else
                            filterMap.put(filterKey, newLine);
                    } else
                        throw new IllegalStateException("Incompatible header type, not VCFFilterHeaderLine: " + line);
                } else if(key.equals("FORMAT")) {
                    if(line instanceof VCFNamedHeaderLine) {
                        String formatKey = ((VCFNamedHeaderLine) line).getName();
                        if(formatMap.containsKey(formatKey))
                            mergeCompoundHeaderLine(formatMap, formatKey, line, logger);
                        else
                            formatMap.put(formatKey, line);
                    }
                } else if(key.equals("SAMPLE")) {
                    String sampleValue = insertCenter(line.getValue(), name);
                    String sampleKey = parseID(line.getValue()) + "." + name;
                    sampleMap.put(sampleKey, new VCFHeaderLine(key, sampleValue));
                } else if(key.equals("vcfProcessLog")) {
                    if(map.containsKey(key)) {
                        //merge new information into vcfProcessLog
                        String currentValue = line.getValue();
                        if(currentValue.contains("InputVCF") &&
                                currentValue.contains("InputVCFSource") &&
                                currentValue.contains("InputVCFVer") &&
                                currentValue.contains("InputVCFParam") &&
                                currentValue.contains("InputVCFgeneAnno")) {
                            VCFHeaderLine oldLine = map.get(key);
                            String oldValue = oldLine.getValue();

                            //original form should be
                            /*
                            ##vcfProcessLog=<InputVCF=<file1.vcf>, 
                            InputVCFSource=<varCaller1>, 
                            InputVCFVer=<1.0>, 
                            InputVCFParam=<a1,c2>,
                            InputVCFgeneAnno=<anno1.gaf>> 
                             */
                            String[] oldValueSplit = oldValue.split(">", 6);
                            String[] currentValueSplit = currentValue.split("<|>", 12);
                            
                            /*
                            oldValueSplit = ["<InputVCF=<file1.vcf", 
                            ",InputVCFSource=<varCaller1", 
                            ",InputVCFVer=<1.0",
                            ",InputVCFParam=<a1,c2",
                            ",InputVCFgeneAnno=<anno1.gaf", ""] 
                              */
                            
                            /*
                            currentValueSplit = ["", "InputVCF=", "file1.vcf", 
                            ",InputVCFSource=", "varCaller1", 
                            ",InputVCFVer=", "1.0",
                            ",InputVCFParam=", "a1,c2",
                            ",InputVCFgeneAnno=", "anno1.gaf", "", ""] 
                             */
                            
                            String newValue = "" +
                                oldValueSplit[0] + "," + currentValueSplit[2] + ">" +
                                oldValueSplit[1] + "," + currentValueSplit[4] + ">" +
                                oldValueSplit[2] + "," + currentValueSplit[6] + ">" +
                                oldValueSplit[3] + "," + currentValueSplit[8] + ">" +
                                oldValueSplit[4] + "," + currentValueSplit[10] + ">" +
                                oldValueSplit[5] + ">";
                            map.put(key, new VCFHeaderLine(key, newValue));
                        } else {
                            throw new IllegalStateException("Incompatible vcfProcessLog," +
                                    "expected InputVCF, InputVCFSource, InputVCFVer, InputVCFParam," +
                                    "and InputVCFgeneAnno but got " + currentValue);
                        }
                    } else {
                        //add information about merging software
                        String oldValue = line.getValue();
                        if(oldValue.contains("InputVCF") &&
                                oldValue.contains("InputVCFSource") &&
                                oldValue.contains("InputVCFVer") &&
                                oldValue.contains("InputVCFParam") &&
                                oldValue.contains("InputVCFgeneAnno")) {
                            map.put(key, line);
                        } else {
                            throw new IllegalStateException("Incompatible vcfProcessLog," +
                                    "expected InputVCF, InputVCFSource, InputVCFVer, InputVCFParam," +
                                    "and InputVCFgeneAnno but got " + oldValue);
                        }
                    }
                } else {
                    if(map.containsKey(key)) {
                        VCFHeaderLine other = map.get(key);
                        if(line.equals(other))
                            continue;
                        else {
                            if ( logger != null ) conflictWarner.warn(line, "Ignoring header line already in map: this header line = " + line + " already present header = " + other);
                            continue;
                        }

                    }else
                        map.put(key, line);

                }

            }
        }
        Set<VCFHeaderLine> headers = new HashSet<VCFHeaderLine>(map.values());
        headers.addAll(infoMap.values());
        headers.addAll(filterMap.values());
        headers.addAll(formatMap.values());
        headers.addAll(sampleMap.values());

        return headers;
    }


    public static String rsIDOfFirstRealVariant(List<VariantContext> VCs, VariantContext.Type type) {
        if ( VCs == null )
            return null;

        String rsID = null;
        for ( VariantContext vc : VCs ) {
            if ( vc.getType() == type ) {
                rsID = vc.getID();
                break;
            }
        }

        return rsID;
    }
}