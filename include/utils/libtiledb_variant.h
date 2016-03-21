/**
 * Copyright (c) 2016 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "query_variants.h"

class Factory {
  private:
    VariantStorageManager *sm;
    VariantQueryProcessor *qp;

    // path to TileDb workspace
    std::string workspace;
    // Name of the array in the workspace
    std::string array_name;

    // Flag to know if the QueryProcessor object needs to be reset
    bool reset_qp;

  public:
    Factory() {
        sm = NULL;
        qp = NULL;
        workspace = "";
        array_name = "";
        reset_qp = true;
    }

    void clear();

    ~Factory() {
        clear();
    }

    //Stats struct
    GTProfileStats stats;

    // Get functions that do lazy initialization of the Tile DB Objects
    VariantStorageManager *getStorageManager(std::string &workspace);
    VariantQueryProcessor *getVariantQueryProcessor(std::string &workspace, 
                                                    const std::string& array_name);
};

extern "C" void db_query_column(std::string workspace, 
                                std::string array_name, 
                                uint64_t query_interval_idx, 
                                Variant& v, 
                                VariantQueryConfig& config); 

extern "C" void db_query_column_range(std::string workspace, 
                                      std::string array_name, 
                                      uint64_t query_interval_idx, 
                                      std::vector<Variant>& variants, 
                                      VariantQueryConfig& query_config, GA4GHPagingInfo* paging_info=0);

extern "C" void db_cleanup();

extern "C" void test_C_pointers(Variant& variant);
