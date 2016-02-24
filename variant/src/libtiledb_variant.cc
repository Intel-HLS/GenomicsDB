#include <iostream>

#include "libtiledb_variant.h"

StorageManager *Factory::getStorageManager(std::string &workspace) {
  if( workspace.compare(this->workspace) != 0 ) {
      if( sm != NULL ) {
          clear();
      }
      // Create storage manager
      // The input is the path to its workspace (the path must exist).
      sm = new StorageManager(workspace);
      this->workspace = workspace;
      // Set reset_qp flag since the sm object has changed
      reset_qp = true;
  }
  return sm;
}

VariantQueryProcessor *Factory::getVariantQueryProcessor(std::string &workspace, const std::string& array_name) {
  if( reset_qp || workspace.compare(this->workspace) != 0 ) {
      // Create query processor
      // The first input is the path to its workspace (the path must exist).
      qp = new VariantQueryProcessor(getStorageManager(workspace), array_name);
      this->workspace = workspace;
      reset_qp = false;
  }
  return qp;
}

void Factory::clear() {
    if( sm == NULL ) {
        return;
    }
    try { 
        sm->close_array(qp->get_array_descriptor());
        delete qp;
    }
    catch (...) { }
    delete sm;
    sm = NULL;
    qp = NULL;
    workspace.clear();
    array_name.clear();
    reset_qp = true;
}

extern "C" void db_query_column(std::string workspace, std::string array_name, 
        uint64_t query_interval_idx, Variant& variant, VariantQueryConfig& query_config) {
    // Init Storage Manager object in the Factory class as 
    // both ArrayDescriptor and Query Processor use it 
    Factory f;
    VariantQueryProcessor *qp = f.getVariantQueryProcessor(workspace, array_name);
    //Do book-keeping, if not already done
    if(!query_config.is_bookkeeping_done())
    {
      qp->do_query_bookkeeping(qp->get_array_schema(), query_config);
      variant = std::move(Variant(&query_config));
      variant.resize_based_on_query();
    }
    f.stats.increment_num_queries();
    qp->gt_get_column(qp->get_array_descriptor(), query_config, query_interval_idx, variant, &f.stats);
}

extern "C" void db_query_column_range(std::string workspace, std::string array_name, 
        uint64_t query_interval_idx, std::vector<Variant>& variants, VariantQueryConfig& query_config, GA4GHPagingInfo* paging_info) {
    // Init Storage Manager object in the Factory class as 
    // both ArrayDescriptor and Query Processor use it 
    Factory f;
    VariantQueryProcessor *qp = f.getVariantQueryProcessor(workspace, array_name);
    //Do book-keeping, if not already done
    if(!query_config.is_bookkeeping_done())
        qp->do_query_bookkeeping(qp->get_array_schema(), query_config);
    qp->gt_get_column_interval(qp->get_array_descriptor(), query_config, query_interval_idx, variants, paging_info, &f.stats);
    if(paging_info == 0 || paging_info->is_query_completed())
        f.stats.increment_num_queries();
}

extern "C" void db_cleanup() {
    // f.clear();
}

template<class T>
void print_vector(T* vec, unsigned size)
{
    std::cout << "[ ";
    for(auto i=0u;i<size;++i)
        std::cout << vec[i] << ",";
    std::cout << " ]\n";
}

#define PRINT_MACRO(X)  print_vector<X>(reinterpret_cast<X*>(ptr), size);

enum TestCPointersEnum
{
  TEST_C_POINTER_INT=0u,
  TEST_C_POINTER_INT64_T,
  TEST_C_POINTER_UNSIGNED,
  TEST_C_POINTER_UINT64_T,
  TEST_C_POINTER_FLOAT,
  TEST_C_POINTER_DOUBLE,
  TEST_C_POINTER_CHAR_PTR,
};

extern "C" void test_C_pointers(Variant& variant)
{
    auto type_to_int = std::unordered_map<std::type_index, unsigned> {
        { std::type_index(typeid(int)), TEST_C_POINTER_INT },
        { std::type_index(typeid(int64_t)), TEST_C_POINTER_INT64_T },
        { std::type_index(typeid(unsigned)), TEST_C_POINTER_UNSIGNED },
        { std::type_index(typeid(uint64_t)), TEST_C_POINTER_UINT64_T },
        { std::type_index(typeid(float)), TEST_C_POINTER_FLOAT },
        { std::type_index(typeid(double)), TEST_C_POINTER_DOUBLE },
        { std::type_index(typeid(char)), TEST_C_POINTER_CHAR_PTR }
    };
    for(Variant::valid_calls_iterator iter=variant.begin();iter!=variant.end();++iter)
    {
        VariantCall& curr_call = *iter;
        for(auto i=0u;i<curr_call.get_num_fields();++i)
        {
            auto& field_ptr = curr_call.get_field(i);   //returns unique_ptr<VariantFieldBase>&
            if(field_ptr.get())
            {
                unsigned size = 0;
                char* ptr = 0;
                bool allocated = false;
                //The function may allocate memory which must be freed by the client code
                //For example, when querying the ALT field, the function will allocate array of N char*
                auto type_index = field_ptr->get_C_pointers(size, reinterpret_cast<void**>(&ptr), allocated);
                if(type_to_int.find(type_index) == type_to_int.end())
                {
                  std::cerr << "Unknown type for field idx "<<i<<", skipping\n";
                  continue;
                }
                switch(type_to_int[type_index])
                {
                    case TEST_C_POINTER_INT:
                        PRINT_MACRO(int);
                        break;
                    case TEST_C_POINTER_INT64_T:
                        PRINT_MACRO(int64_t);
                        break;  
                    case TEST_C_POINTER_UNSIGNED:
                        PRINT_MACRO(unsigned);
                        break;
                    case TEST_C_POINTER_UINT64_T:
                        PRINT_MACRO(uint64_t);
                        break;
                    case TEST_C_POINTER_FLOAT:
                        PRINT_MACRO(float);
                        break;
                    case TEST_C_POINTER_DOUBLE:
                        PRINT_MACRO(double);
                        break;
                    case TEST_C_POINTER_CHAR_PTR:
                        PRINT_MACRO(char*);
                        break;
                }
                if(allocated)
                {
                  if(size > 1)
                    delete[] ptr;
                  else
                    delete ptr;
                }
            }
        }
    }
}
