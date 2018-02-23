/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2018 Intel Corporation
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package com.intel.genomicsdb.reader.model;

import org.json.simple.JSONObject;

/**
 * Represents parameters to use in a GenomicsDB query.
 */
public class QueryParams {
    private static final boolean SCAN_FULL = true;
    private String workspace;
    private String arrayName;
    private String vidMappingFile;
    private String callsetMappingFile;
    private String referenceGenome;
    private String templateVCFHeaderFilename;
    private Boolean produceGtField;

    public QueryParams(final String workspace, final String referenceGenome, final String templateVCFHeaderFilename) {
        this(workspace, "", null, null, referenceGenome, templateVCFHeaderFilename, false);
    }

    public QueryParams(final String workspace, final String vidMappingFile, final String callsetMappingFile,
                       final String referenceGenome, final String templateVCFHeaderFilename, final Boolean produceGtField) {
        this(workspace, "", vidMappingFile, callsetMappingFile, referenceGenome, templateVCFHeaderFilename, produceGtField);
    }

    public QueryParams(final QueryParams source) {
        this(source.workspace, source.arrayName, source.vidMappingFile, source.callsetMappingFile,
                source.referenceGenome, source.templateVCFHeaderFilename, source.produceGtField);
    }

    private QueryParams(final String workspace, final String arrayName, final String vidMappingFile,
                       final String callsetMappingFile, final String referenceGenome,
                       final String templateVCFHeaderFilename, final Boolean produceGtField) {
        this.workspace = workspace;
        this.arrayName = arrayName;
        this.vidMappingFile = vidMappingFile;
        this.callsetMappingFile = callsetMappingFile;
        this.referenceGenome = referenceGenome;
        this.templateVCFHeaderFilename = templateVCFHeaderFilename;
        this.produceGtField = produceGtField;
    }

    public String getWorkspace() {
        return workspace;
    }

    public String getArrayName() {
        return arrayName;
    }

    public void setArrayName(String arrayName) {
        this.arrayName = arrayName;
    }

    /**
     * Create a JSON string based on the attributes of this class.
     * @return A JSON string
     */
    public String toJsonString(){
        JSONObject queryParams = new JSONObject();
        queryParams.put("scan_full", SCAN_FULL);
        queryParams.put("workspace", this.workspace);
        queryParams.put("array", this.arrayName);
        if (this.vidMappingFile != null) queryParams.put("vid_mapping_file", this.vidMappingFile);
        if (this.callsetMappingFile != null) queryParams.put("callset_mapping_file", this.callsetMappingFile);
        queryParams.put("reference_genome", this.referenceGenome);
        if (this.templateVCFHeaderFilename != null) queryParams.put("templateVCFHeaderFilename", this.templateVCFHeaderFilename);
        queryParams.put("produceGtField", this.produceGtField);
        return queryParams.toJSONString();
    }
}
