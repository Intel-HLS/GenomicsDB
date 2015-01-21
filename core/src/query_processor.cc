/**
 * @file   query_processor.cc
 * @author Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * @section LICENSE
 *
 * The MIT License
 * 
 * @copyright Copyright (c) 2014 Stavros Papadopoulos <stavrosp@csail.mit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @section DESCRIPTION
 *
 * This file implements the QueryProcessor class.
 */
  
#include "query_processor.h"
#include <stdio.h>
#include <typeinfo>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>

/* ------------------- Genotyping functions ------------------ */

/******************************************************
********************* CONSTRUCTORS ********************
******************************************************/

QueryProcessor::GTColumn::GTColumn(int64_t col, uint64_t row_num) {
  ALT_.resize(row_num);
  col_ = col;
  REF_.resize(row_num);
  PL_.resize(row_num);
}

/******************************************************
*********************** OPERATIONS ********************
******************************************************/

// TODO: Implement function for deriving the final ALT and PL values

// TODO: Implement the genotyping function

/* ----------------- QueryProcessor functions ---------------- */

/******************************************************
********************* CONSTRUCTORS ********************
******************************************************/

QueryProcessor::QueryProcessor(const std::string& workspace, 
                               StorageManager& storage_manager) 
    : storage_manager_(storage_manager) {
  set_workspace(workspace);
  create_workspace(); 
}

void QueryProcessor::export_to_CSV(const StorageManager::ArrayDescriptor* ad,
                                   const std::string& filename) const { 
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::string& array_name = array_schema.array_name();
  const unsigned int attribute_num = array_schema.attribute_num();
  const unsigned int dim_num = array_schema.dim_num();
  
  // Prepare CSV file
  CSVFile csv_file(filename, CSVFile::WRITE);

  // Create and initialize tile iterators
  StorageManager::const_iterator *tile_its = 
      new StorageManager::const_iterator[attribute_num+1];
  StorageManager::const_iterator tile_it_end;
  initialize_tile_its(ad, tile_its, tile_it_end);

  // Create cell iterators
  Tile::const_iterator* cell_its = 
      new Tile::const_iterator[attribute_num+1];
  Tile::const_iterator cell_it_end;

  // Iterate over all tiles
  while(tile_its[attribute_num] != tile_it_end) {
    // Iterate over all cells of each tile
    initialize_cell_its(tile_its, attribute_num, cell_its, cell_it_end);

    while(cell_its[attribute_num] != cell_it_end) { 
      csv_file << cell_to_csv_line(cell_its, attribute_num);
      advance_cell_its(attribute_num, cell_its);
    }
 
    advance_tile_its(attribute_num, tile_its);
  }

  // Clean up 
  delete [] tile_its;
  delete [] cell_its;
}

QueryProcessor::GTColumn* QueryProcessor::gt_get_column(
    const StorageManager::ArrayDescriptor* ad, uint64_t col) const {
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  const std::vector<std::pair<double, double> >& dim_domains =
      array_schema.dim_domains();
  uint64_t row_num = dim_domains[0].second - dim_domains[0].first + 1;
  unsigned int attribute_num = array_schema.attribute_num();

  // Check that column falls into the domain of the second dimension
  assert(col >= dim_domains[1].first && col <= dim_domains[1].second);

  // Indicates how many rows have been filled.
  uint64_t filled_rows = 0;

  // Initialize reverse tile iterators for 
  // END, REF, ALT, PL, NULL, OFFSETS, coordinates
  // The reverse tile iterator will start with the tiles
  // of the various attributes that have the largest
  // id that either intersect with col, or precede col.
  StorageManager::const_reverse_iterator* tile_its;
  StorageManager::const_reverse_iterator tile_it_end;  
  unsigned int gt_attribute_num = 
      gt_initialize_tile_its(ad, tile_its, tile_it_end, col);

  // Create and initialize GenotypingColumn members
  QueryProcessor::GTColumn* gt_column = new GTColumn(col, row_num);

  // Create cell iterators
  Tile::const_reverse_iterator cell_it, cell_it_end;

  // Fill the genotyping column
  while(tile_its[gt_attribute_num] != tile_it_end && filled_rows < row_num) {
    // Initialize cell iterators for the coordinates
    cell_it = (*tile_its[gt_attribute_num]).rbegin();
    cell_it_end = (*tile_its[gt_attribute_num]).rend();
    while(cell_it != cell_it_end && filled_rows < row_num) {
      std::vector<int64_t> next_coord = *cell_it;
      // If next cell is not on the right of col, and corresponds to 
      // uninvestigated row
      if(next_coord[1] <= col && gt_column->REF_[next_coord[0]] == "") {
        gt_fill_row(gt_column, next_coord[0], cell_it.pos(), tile_its);
        ++filled_rows;
      }
      ++cell_it;
    }
    advance_tile_its(gt_attribute_num, tile_its);
  }

  assert(filled_rows == row_num);

  delete [] tile_its;

  return gt_column;
}

void QueryProcessor::join(const StorageManager::ArrayDescriptor* ad_A, 
                          const StorageManager::ArrayDescriptor* ad_B,
                          const std::string& result_array_name) const {
  // For easy reference
  const ArraySchema& array_schema_A = ad_A->array_info()->array_schema_;
  const ArraySchema& array_schema_B = ad_B->array_info()->array_schema_;

  std::pair<bool,std::string> can_join = 
      ArraySchema::join_compatible(array_schema_A, array_schema_B);

  if(!can_join.first)
    throw QueryProcessorException(std::string("[QueryProcessor] Input arrays "
                                  " are not join-compatible.") + 
                                  can_join.second);

  ArraySchema array_schema_C = ArraySchema::create_join_result_schema(
                                   array_schema_A, 
                                   array_schema_B, 
                                   result_array_name);
  
  if(array_schema_A.has_regular_tiles())
    join_regular(ad_A, ad_B, array_schema_C);
  else 
    join_irregular(ad_A, ad_B, array_schema_C);
} 

void QueryProcessor::subarray(const StorageManager::ArrayDescriptor* ad,
                              const Tile::Range& range,
                              const std::string& result_array_name) const { 
  if(ad->array_schema().has_regular_tiles())
    subarray_regular(ad, range, result_array_name);
  else 
    subarray_irregular(ad, range, result_array_name);
}

/******************************************************
******************* PRIVATE METHODS *******************
******************************************************/

inline
void QueryProcessor::advance_cell_its(unsigned int attribute_num,
                                      Tile::const_iterator* cell_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
      ++cell_its[i];
}

inline
void QueryProcessor::advance_cell_its(unsigned int attribute_num,
                                      Tile::const_iterator* cell_its,
                                      int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    cell_its[i] += step;
}

inline
void QueryProcessor::advance_cell_its(
    unsigned int attribute_num, Tile::const_reverse_iterator* cell_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
      ++cell_its[i];
}

inline
void QueryProcessor::advance_cell_its(unsigned int attribute_num,
                                      Tile::const_reverse_iterator* cell_its,
                                      int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    cell_its[i] += step;
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_iterator* tile_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
    ++tile_its[i];
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_iterator* tile_its, 
    int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    tile_its[i] += step;
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_reverse_iterator* tile_its) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
    ++tile_its[i];
}

inline
void QueryProcessor::advance_tile_its(
    unsigned int attribute_num, 
    StorageManager::const_reverse_iterator* tile_its, 
    int64_t step) const {
  for(unsigned int i=0; i<attribute_num; i++) 
    tile_its[i] += step;
}

inline
void QueryProcessor::append_cell(const Tile::const_iterator* cell_its,
                                 Tile** tiles,
                                 unsigned int attribute_num) const {
  for(unsigned int i=0; i<=attribute_num; i++) 
    *tiles[i] << cell_its[i]; 
}

inline
void QueryProcessor::append_cell(const Tile::const_iterator* cell_its_A,
                                 const Tile::const_iterator* cell_its_B,
                                 Tile** tiles_C,
                                 unsigned int attribute_num_A,
                                 unsigned int attribute_num_B) const {
  for(unsigned int i=0; i<attribute_num_A; i++) 
    *tiles_C[i] << cell_its_A[i]; 
  for(unsigned int i=0; i<=attribute_num_B; i++)
    *tiles_C[attribute_num_A+i] << cell_its_B[i]; 
}

inline
CSVLine QueryProcessor::cell_to_csv_line(const Tile::const_iterator* cell_its,
                                         unsigned int attribute_num) const {
  CSVLine csv_line;

  // Append coordinates first
  cell_its[attribute_num] >> csv_line;
  // Append attribute values next
  for(unsigned int i=0; i<attribute_num; i++)
    cell_its[i] >> csv_line;

  return csv_line;
}

void QueryProcessor::create_workspace() const {
  struct stat st;
  stat(workspace_.c_str(), &st);

  // If the workspace does not exist, create it
  if(!S_ISDIR(st.st_mode)) { 
    int dir_flag = mkdir(workspace_.c_str(), S_IRWXU);
    assert(dir_flag == 0);
  }
}

inline
void QueryProcessor::get_tiles(
    const StorageManager::ArrayDescriptor* ad, 
    uint64_t tile_id, const Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Get attribute tiles
  for(unsigned int i=0; i<=attribute_num; i++) 
   tiles[i] = storage_manager_.get_tile(ad, i, tile_id);
}

bool QueryProcessor::path_exists(const std::string& path) const {
  struct stat st;
  stat(path.c_str(), &st);
  return S_ISDIR(st.st_mode);
}

void QueryProcessor::gt_fill_row(
    GTColumn* gt_column, int64_t row, int64_t pos,
    const StorageManager::const_reverse_iterator* tile_its) const {
  // First check if the row is NULL
  int64_t END_v = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[0]).cell(pos);
  if(END_v < gt_column->col_) {
    gt_column->REF_[row] = "$";
    return;
  }

  // Retrieve the offsets
  const AttributeTile<int64_t>& OFFSETS_tile = 
      static_cast<const AttributeTile<int64_t>& >(*tile_its[5]);
  int64_t REF_offset = OFFSETS_tile.cell(pos*5);
  int64_t ALT_offset = OFFSETS_tile.cell(pos*5+1);
  int64_t PL_offset = OFFSETS_tile.cell(pos*5+4);

  // Retrieve the NULL bitmap
  const AttributeTile<int>& NULL_tile = 
      static_cast<const AttributeTile<int>& >(*tile_its[4]);
  int NULL_bitmap = NULL_tile.cell(pos);

  char c;
  int i;

  // Fill the REF
  const AttributeTile<char>& REF_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[1]);
  std::string REF_s = "";
  i = 0;
  while((c = REF_tile.cell(REF_offset+i)) != '\0') { 
    REF_s.push_back(c);
    ++i;
  }
  gt_column->REF_[row] = REF_s;

  // Fill the ALT values
  const AttributeTile<char>& ALT_tile = 
      static_cast<const AttributeTile<char>& >(*tile_its[2]);
  i = 0;
  std::string ALT_s = "";
  while((c = ALT_tile.cell(ALT_offset+i)) != '&') {
    if(c == '\0') {
      gt_column->ALT_[row].push_back(ALT_s);
      ALT_s = "";
    } else {
      ALT_s.push_back(c);
    }
    i++;
  }
  assert(ALT_s == "");
  gt_column->ALT_[row].push_back("&");
   
  // Fill the PL values
  if(NULL_bitmap & 1 == 0) { // If the PL values are not NULL
    const AttributeTile<int>& PL_tile = 
        static_cast<const AttributeTile<int>& >(*tile_its[3]);
    int ALT_num = gt_column->ALT_[row].size(); 
    int PL_num = (ALT_num+1)*(ALT_num+2)/2;
    for(int i=0; i<PL_num; i++) 
      gt_column->PL_[row].push_back(PL_tile.cell(PL_offset+i));
  }
}

inline
unsigned int QueryProcessor::gt_initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    StorageManager::const_reverse_iterator*& tile_its, 
    StorageManager::const_reverse_iterator& tile_it_end,
    uint64_t col) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Create reverse iterators
  tile_its = new StorageManager::const_reverse_iterator[7];
  // Find the rank of the tile the left sweep starts from.
  uint64_t start_rank = storage_manager_.get_left_sweep_start_rank(ad, col);

  // END
  tile_its[0] = storage_manager_.rbegin(ad, 0, start_rank);
  // REF
  tile_its[1] = storage_manager_.rbegin(ad, 1, start_rank);
  // ALT
  tile_its[2] = storage_manager_.rbegin(ad, 2, start_rank);
  // PL
  tile_its[3] = storage_manager_.rbegin(ad, 20, start_rank);
  // NULL
  tile_its[4] = storage_manager_.rbegin(ad, 21, start_rank);
  // OFFSETS
  tile_its[5] = storage_manager_.rbegin(ad, 22, start_rank);
  // coordinates
  tile_its[6] = storage_manager_.rbegin(ad, attribute_num, start_rank);
  tile_it_end = storage_manager_.rend(ad, attribute_num);

  // The number of attributes is 6, and the coordinates is the extra one
  return 6;
}

inline
void QueryProcessor::initialize_cell_its(
    const Tile** tiles, unsigned int attribute_num,
    Tile::const_iterator* cell_its, Tile::const_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = tiles[i]->begin();
  cell_it_end = tiles[attribute_num]->end();
}

inline
void QueryProcessor::initialize_cell_its(
    const StorageManager::const_iterator* tile_its, unsigned int attribute_num,
    Tile::const_iterator* cell_its, Tile::const_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = (*tile_its[i]).begin();
  cell_it_end = (*tile_its[attribute_num]).end();
}

inline
void QueryProcessor::initialize_cell_its(
    const StorageManager::const_reverse_iterator* tile_its, 
    unsigned int attribute_num,
    Tile::const_iterator* cell_its, 
    Tile::const_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = (*tile_its[i]).begin();
  cell_it_end = (*tile_its[attribute_num]).end();
}

inline
void QueryProcessor::initialize_cell_its(
    const StorageManager::const_reverse_iterator* tile_its, 
    unsigned int attribute_num,
    Tile::const_reverse_iterator* cell_its, 
    Tile::const_reverse_iterator& cell_it_end) const {
  for(unsigned int i=0; i<=attribute_num; i++)
    cell_its[i] = (*tile_its[i]).rbegin();
  cell_it_end = (*tile_its[attribute_num]).rend();
}

inline
void QueryProcessor::initialize_tile_its(
    const StorageManager::ArrayDescriptor* ad,
    StorageManager::const_iterator* tile_its, 
    StorageManager::const_iterator& tile_it_end) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  for(unsigned int i=0; i<=attribute_num; i++)
    tile_its[i] = storage_manager_.begin(ad, i);
  tile_it_end = storage_manager_.end(ad, attribute_num);
}

void QueryProcessor::join_irregular(const StorageManager::ArrayDescriptor* ad_A, 
                                    const StorageManager::ArrayDescriptor* ad_B,
                                    const ArraySchema& array_schema_C) const {
  // For easy reference
  const ArraySchema& array_schema_A = ad_A->array_schema();
  const ArraySchema& array_schema_B = ad_B->array_schema();
  unsigned int attribute_num_A = array_schema_A.attribute_num();
  unsigned int attribute_num_B = array_schema_B.attribute_num();
  unsigned int attribute_num_C = array_schema_C.attribute_num();

  // Prepare result array
  const StorageManager::ArrayDescriptor* ad_C = 
      storage_manager_.open_array(array_schema_C);

  // Create tiles 
  const Tile** tiles_A = new const Tile*[attribute_num_A+1];
  const Tile** tiles_B = new const Tile*[attribute_num_B+1];
  Tile** tiles_C = new Tile*[attribute_num_C+1];
  
  // Create and initialize tile iterators
  StorageManager::const_iterator *tile_its_A = 
      new StorageManager::const_iterator[attribute_num_A+1];
  StorageManager::const_iterator *tile_its_B = 
      new StorageManager::const_iterator[attribute_num_B+1];
  StorageManager::const_iterator tile_it_end_A;
  StorageManager::const_iterator tile_it_end_B;
  initialize_tile_its(ad_A, tile_its_A, tile_it_end_A);
  initialize_tile_its(ad_B, tile_its_B, tile_it_end_B);
  
  // Create cell iterators
  Tile::const_iterator* cell_its_A = 
      new Tile::const_iterator[attribute_num_A+1];
  Tile::const_iterator* cell_its_B = 
      new Tile::const_iterator[attribute_num_B+1];
  Tile::const_iterator cell_it_end_A, cell_it_end_B;

  // Auxiliary variables storing the number of skipped tiles when joining.
  // It is used to advance only the coordinates iterator when a tile is
  // finished/skipped, and then efficiently advance the attribute iterators only 
  // when a tile joins.
  int64_t skipped_tiles_A = 0;
  int64_t skipped_tiles_B = 0;

  // Initialize tiles with id 0 for C (result array)
  new_tiles(array_schema_C, 0, tiles_C); 

  // Join algorithm
  while(tile_its_A[attribute_num_A] != tile_it_end_A &&
        tile_its_B[attribute_num_B] != tile_it_end_B) {
    // Potential join result generation
    if(may_join(tile_its_A[attribute_num_A], tile_its_B[attribute_num_B])) {
      // Update iterators in A
      if(skipped_tiles_A) {
        advance_tile_its(attribute_num_A, tile_its_A, skipped_tiles_A);
        skipped_tiles_A = 0;
        initialize_cell_its(tile_its_A, attribute_num_A, 
                            cell_its_A, cell_it_end_A);
      }
      // Update iterators in B
      if(skipped_tiles_B) {
        advance_tile_its(attribute_num_B, tile_its_B, skipped_tiles_B);
        skipped_tiles_B = 0;
        initialize_cell_its(tile_its_B, attribute_num_B, 
                            cell_its_B, cell_it_end_B);
      }
      // Join the tiles
      join_tiles_irregular(attribute_num_A, cell_its_A, cell_it_end_A, 
                           attribute_num_B, cell_its_B, cell_it_end_B,
                           ad_C, tiles_C);
    }

    // Check which tile precedes the other in the global order
    // Note that operator '<', when the operands are from different
    // arrays, returns true if the first tile precedes the second
    // in the global order by checking their bounding coordinates.
    if(tile_its_A[attribute_num_A] < tile_its_B[attribute_num_B]) {
      ++tile_its_A[attribute_num_A];
      ++skipped_tiles_A;
    }
    else {
      ++tile_its_B[attribute_num_B];
      ++skipped_tiles_B;
    }
  }
  
  // Send the lastly created tiles to storage manager
  store_tiles(ad_C, tiles_C);

  // Close result array
  storage_manager_.close_array(ad_C);

  // Clean up
  delete [] tiles_A;
  delete [] tiles_B;
  delete [] tiles_C;
  delete [] tile_its_A;
  delete [] tile_its_B;
  delete [] cell_its_A;
  delete [] cell_its_B;
}

void QueryProcessor::join_regular(const StorageManager::ArrayDescriptor* ad_A, 
                                  const StorageManager::ArrayDescriptor* ad_B,
                                  const ArraySchema& array_schema_C) const {
  // For easy reference
  const ArraySchema& array_schema_A = ad_A->array_schema();
  const ArraySchema& array_schema_B = ad_B->array_schema();
  unsigned int attribute_num_A = array_schema_A.attribute_num();
  unsigned int attribute_num_B = array_schema_B.attribute_num();
  unsigned int attribute_num_C = array_schema_C.attribute_num();
  uint64_t tile_id_A, tile_id_B;

  // Prepare result array
  const StorageManager::ArrayDescriptor* ad_C = 
      storage_manager_.open_array(array_schema_C);

  // Create tiles 
  const Tile** tiles_A = new const Tile*[attribute_num_A+1];
  const Tile** tiles_B = new const Tile*[attribute_num_B+1];
  Tile** tiles_C = new Tile*[attribute_num_C+1];
  
  // Create and initialize tile iterators
  StorageManager::const_iterator *tile_its_A = 
      new StorageManager::const_iterator[attribute_num_A+1];
  StorageManager::const_iterator *tile_its_B = 
      new StorageManager::const_iterator[attribute_num_B+1];
  StorageManager::const_iterator tile_it_end_A;
  StorageManager::const_iterator tile_it_end_B;
  initialize_tile_its(ad_A, tile_its_A, tile_it_end_A);
  initialize_tile_its(ad_B, tile_its_B, tile_it_end_B);
  
  // Create cell iterators
  Tile::const_iterator* cell_its_A = 
      new Tile::const_iterator[attribute_num_A+1];
  Tile::const_iterator* cell_its_B = 
      new Tile::const_iterator[attribute_num_B+1];
  Tile::const_iterator cell_it_end_A, cell_it_end_B;

  // Auxiliary variables storing the number of skipped tiles when joining.
  // It is used to advance only the coordinates iterator when a tile is
  // finished/skipped, and then efficiently advance the attribute iterators only
  // when a tile joins.
  int64_t skipped_tiles_A = 0;
  int64_t skipped_tiles_B = 0;

  // Join algorithm
  while(tile_its_A[attribute_num_A] != tile_it_end_A &&
        tile_its_B[attribute_num_B] != tile_it_end_B) {
    tile_id_A = tile_its_A[attribute_num_A].tile_id();
    tile_id_B = tile_its_B[attribute_num_B].tile_id();

    // Potential join result generation
    if(tile_id_A == tile_id_B) {
      // Update iterators in A
      if(skipped_tiles_A) {
        advance_tile_its(attribute_num_A, tile_its_A, skipped_tiles_A);
        skipped_tiles_A = 0;
        initialize_cell_its(tile_its_A, attribute_num_A, 
                            cell_its_A, cell_it_end_A);
      }
      // Update iterators in B
      if(skipped_tiles_B) {
        advance_tile_its(attribute_num_B, tile_its_B, skipped_tiles_B);
        skipped_tiles_B = 0;
        initialize_cell_its(tile_its_B, attribute_num_B, 
                            cell_its_B, cell_it_end_B);
      }

      // Initialize tiles for C (result array)
      new_tiles(array_schema_C, tile_id_A, tiles_C);
      // Join the tiles
      join_tiles_regular(attribute_num_A, cell_its_A, cell_it_end_A, 
                         attribute_num_B, cell_its_B, cell_it_end_B,
                         ad_C, tiles_C);
      // Send the created tiles to storage manager
      store_tiles(ad_C, tiles_C);
    }

    // Tile precedence in the case of regular tiles is simply determined
    // by the order of the tile ids.
    if(tile_id_A < tile_id_B) {
      ++tile_its_A[attribute_num_A];
      ++skipped_tiles_A;
    } else if(tile_id_A > tile_id_B) {
      ++tile_its_B[attribute_num_B];
      ++skipped_tiles_B;
    } else { // tile_id_A == tile_id_B, advance both
      ++tile_its_A[attribute_num_A];
      ++skipped_tiles_A;
      ++tile_its_B[attribute_num_B];
      ++skipped_tiles_B;
    }
  }
  
  // Close result array
  storage_manager_.close_array(ad_C);

  // Clean up
  delete [] tiles_A;
  delete [] tiles_B;
  delete [] tiles_C;
  delete [] tile_its_A;
  delete [] tile_its_B;
  delete [] cell_its_A;
  delete [] cell_its_B;
}

void QueryProcessor::join_tiles_irregular(
    unsigned int attribute_num_A, Tile::const_iterator* cell_its_A,
    Tile::const_iterator& cell_it_end_A, 
    unsigned int attribute_num_B, Tile::const_iterator* cell_its_B,
    Tile::const_iterator& cell_it_end_B,
    const StorageManager::ArrayDescriptor* ad_C, Tile** tiles_C) const {
  // For easy reference
  const ArraySchema& array_schema_C = ad_C->array_schema();
  uint64_t capacity = array_schema_C.capacity();
  unsigned int attribute_num_C = array_schema_C.attribute_num();

  // Auxiliary variables storing the number of skipped cells when joining.
  // It is used to advance only the coordinates iterator when a cell is
  // finished/skipped, and then efficiently advance the attribute iterators only 
  // when a cell joins.
  int64_t skipped_cells_A = 0;
  int64_t skipped_cells_B = 0;

  while(cell_its_A[attribute_num_A] != cell_it_end_A &&
        cell_its_B[attribute_num_B] != cell_it_end_B) {
    // If the coordinates are equal
    // Note that operator '==', when the operands correspond to different
    // tiles, returns true if the cell values pointed by the iterators
    // are equal.
    if(cell_its_A[attribute_num_A] == cell_its_B[attribute_num_B]) {      
      advance_cell_its(attribute_num_A, cell_its_A, skipped_cells_A);
      advance_cell_its(attribute_num_B, cell_its_B, skipped_cells_B);
      skipped_cells_A = 0;
      skipped_cells_B = 0;
      if(tiles_C[attribute_num_C]->cell_num() == capacity) {
        uint64_t new_tile_id = tiles_C[attribute_num_C]->tile_id() + 1;
        store_tiles(ad_C, tiles_C);
        new_tiles(array_schema_C, new_tile_id, tiles_C); 
      }
      append_cell(cell_its_A, cell_its_B, tiles_C, 
                  attribute_num_A, attribute_num_B);
      advance_cell_its(attribute_num_A, cell_its_A);
      advance_cell_its(attribute_num_B, cell_its_B);
    // Otherwise check which cell iterator to advance
    } else {
      if(array_schema_C.precedes(cell_its_A[attribute_num_A],
                                 cell_its_B[attribute_num_B])) {
        ++cell_its_A[attribute_num_A];
        ++skipped_cells_A;
      } else {
        ++cell_its_B[attribute_num_B];
        ++skipped_cells_B;
      }
    }
  }
}

void QueryProcessor::join_tiles_regular(
    unsigned int attribute_num_A, Tile::const_iterator* cell_its_A,
    Tile::const_iterator& cell_it_end_A, 
    unsigned int attribute_num_B, Tile::const_iterator* cell_its_B,
    Tile::const_iterator& cell_it_end_B,
    const StorageManager::ArrayDescriptor* ad_C, Tile** tiles_C) const {
  // For easy reference
  const ArraySchema& array_schema_C = ad_C->array_schema();
  unsigned int attribute_num_C = array_schema_C.attribute_num();

  // Auxiliary variables storing the number of skipped cells when joining.
  // It is used to advance only the coordinates iterator when a cell is
  // finished/skipped, and then efficiently advance the attribute iterators only 
  // when a cell joins.
  int64_t skipped_cells_A = 0;
  int64_t skipped_cells_B = 0;

  while(cell_its_A[attribute_num_A] != cell_it_end_A &&
        cell_its_B[attribute_num_B] != cell_it_end_B) {
    // If the coordinates are equal
    // Note that operator '==', when the operands correspond to different
    // tiles, returns true if the cell values pointed by the iterators
    // are equal.
    if(cell_its_A[attribute_num_A] == cell_its_B[attribute_num_B]) {      
      advance_cell_its(attribute_num_A, cell_its_A, skipped_cells_A);
      advance_cell_its(attribute_num_B, cell_its_B, skipped_cells_B);
      skipped_cells_A = 0;
      skipped_cells_B = 0;
      append_cell(cell_its_A, cell_its_B, tiles_C, 
                  attribute_num_A, attribute_num_B);
      advance_cell_its(attribute_num_A, cell_its_A);
      advance_cell_its(attribute_num_B, cell_its_B);
    // Otherwise check which cell iterator to advance
    } else {
      if(array_schema_C.precedes(cell_its_A[attribute_num_A],
                                 cell_its_B[attribute_num_B])) {
        ++cell_its_A[attribute_num_A];
        ++skipped_cells_A;
      } else {
        ++cell_its_B[attribute_num_B];
        ++skipped_cells_B;
      }
    }
  }
}

bool QueryProcessor::may_join(
    const StorageManager::const_iterator& it_A,
    const StorageManager::const_iterator& it_B) const {
  // For easy reference
  const ArraySchema& array_schema_A = it_A.array_schema();
  const MBR& mbr_A = it_A.mbr();
  const MBR& mbr_B = it_B.mbr();

  // Check if the tile MBRs overlap
  if(array_schema_A.has_irregular_tiles() && !overlap(mbr_A, mbr_B))
    return false;

  // For easy reference 
  BoundingCoordinatesPair bounding_coordinates_A = it_A.bounding_coordinates();
  BoundingCoordinatesPair bounding_coordinates_B = it_B.bounding_coordinates();

  // Check if the cell id ranges (along the global order) intersect
  if(!overlap(bounding_coordinates_A, bounding_coordinates_B, array_schema_A))
    return false;

  return true;
}

inline
void QueryProcessor::new_tiles(const ArraySchema& array_schema, 
                               uint64_t tile_id, Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = array_schema.attribute_num();
  uint64_t capacity = array_schema.capacity();

  for(unsigned int i=0; i<=attribute_num; i++)
    tiles[i] = storage_manager_.new_tile(array_schema, i, tile_id, capacity);
}

bool QueryProcessor::overlap(const MBR& mbr_A, const MBR& mbr_B) const {
  assert(mbr_A.size() == mbr_B.size());
  assert(mbr_A.size() % 2 == 0);

  // For easy rederence
  unsigned int dim_num = mbr_A.size() / 2;

  for(unsigned int i=0; i<dim_num; i++) 
    if(mbr_A[2*i+1] < mbr_B[2*i] || mbr_A[2*i] > mbr_B[2*i+1])
      return false;

  return true;
}

bool QueryProcessor::overlap(
    const BoundingCoordinatesPair& bounding_coordinates_A,
    const BoundingCoordinatesPair& bounding_coordinates_B, 
    const ArraySchema& array_schema) const {
  if(array_schema.precedes(bounding_coordinates_A.second, 
                           bounding_coordinates_B.first) ||
     array_schema.succeeds(bounding_coordinates_A.first, 
                           bounding_coordinates_B.second))
    return false;
  else
    return true;
}

inline
void QueryProcessor::set_workspace(const std::string& path) {
  workspace_ = path;
  
  // Replace '~' with the absolute path
  if(workspace_[0] == '~') {
    workspace_ = std::string(getenv("HOME")) +
                 workspace_.substr(1, workspace_.size()-1);
  }

  // Check if the input path is an existing directory 
  assert(path_exists(workspace_));
 
  workspace_ += "/QueryProcessor";
}

inline
void QueryProcessor::store_tiles(const StorageManager::ArrayDescriptor* ad,
                                 Tile** tiles) const {
  // For easy reference
  unsigned int attribute_num = ad->array_schema().attribute_num();

  // Append attribute tiles
  for(unsigned int i=0; i<=attribute_num; i++)
    storage_manager_.append_tile(tiles[i], ad, i);
} 

void QueryProcessor::subarray_irregular(
    const StorageManager::ArrayDescriptor* ad,
    const Tile::Range& range, const std::string& result_array_name) const { 
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  unsigned int attribute_num = array_schema.attribute_num();
  uint64_t capacity = array_schema.capacity();

  // Create tiles
  const Tile** tiles = new const Tile*[attribute_num+1];
  Tile** result_tiles = new Tile*[attribute_num+1];

  // Create cell iterators
  Tile::const_iterator* cell_its = 
      new Tile::const_iterator[attribute_num+1];
  Tile::const_iterator cell_it_end;

  // Prepare result array
  ArraySchema result_array_schema = array_schema.clone(result_array_name);
  const StorageManager::ArrayDescriptor* result_ad = 
      storage_manager_.open_array(result_array_schema);
  
  // Get the tile ids that overlap with the range
  std::vector<std::pair<uint64_t, bool> > overlapping_tile_ids;
  storage_manager_.get_overlapping_tile_ids(ad, range, &overlapping_tile_ids);

  // Initialize tile iterators
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it = 
      overlapping_tile_ids.begin();
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it_end =
      overlapping_tile_ids.end();
      
  // Create result tiles and load input array tiles 
  uint64_t tile_id = 0;
  new_tiles(result_array_schema, tile_id, result_tiles); 

  // Auxiliary variable storing the number of skipped cells when investigating
  // a tile partially overlapping the range. It is used to advance only the
  // coordinates iterator when a cell is not in range, and then efficiently
  // advance the attribute iterators only when a cell falls in the range.
  int64_t skipped;

  // Iterate over all tiles
  for(; tile_id_it != tile_id_it_end; ++tile_id_it) {
    get_tiles(ad, tile_id_it->first, tiles);
    initialize_cell_its(tiles, attribute_num, cell_its, cell_it_end); 
    skipped = 0;

    if(tile_id_it->second) { // Full overlap
      while(cell_its[attribute_num] != cell_it_end) {
        if(result_tiles[attribute_num]->cell_num() == capacity) {
          store_tiles(result_ad, result_tiles);
          new_tiles(result_array_schema, ++tile_id, result_tiles); 
        }
        append_cell(cell_its, result_tiles, attribute_num);
        advance_cell_its(attribute_num, cell_its);
      }
    } else { // Partial overlap
      while(cell_its[attribute_num] != cell_it_end) {
        if(cell_its[attribute_num].cell_inside_range(range)) {
          if(result_tiles[attribute_num]->cell_num() == capacity) {
            store_tiles(result_ad, result_tiles);
            new_tiles(result_array_schema, ++tile_id, result_tiles); 
          }
          advance_cell_its(attribute_num, cell_its, skipped);
          skipped = 0;
          append_cell(cell_its, result_tiles, attribute_num);
          advance_cell_its(attribute_num, cell_its);
        } else { // Advance only the coordinates cell iterator
          skipped++;
          ++cell_its[attribute_num];
        }
      }
    }
  } 

  // Send the lastly created tiles to storage manager
  store_tiles(result_ad, result_tiles);
  
  // Close result array
  storage_manager_.close_array(result_ad);

  // Clean up 
  delete [] tiles;
  delete [] result_tiles;
  delete [] cell_its;
}

void QueryProcessor::subarray_regular(
    const StorageManager::ArrayDescriptor* ad,
    const Tile::Range& range, const std::string& result_array_name) const { 
  // For easy reference
  const ArraySchema& array_schema = ad->array_schema();
  unsigned int attribute_num = array_schema.attribute_num();

  // Create tiles 
  const Tile** tiles = new const Tile*[attribute_num+1];
  Tile** result_tiles = new Tile*[attribute_num+1];

  // Create cell iterators
  Tile::const_iterator* cell_its = 
      new Tile::const_iterator[attribute_num+1];
  Tile::const_iterator cell_it_end;

  // Prepare result array
  ArraySchema result_array_schema = array_schema.clone(result_array_name);
  const StorageManager::ArrayDescriptor* result_ad = 
      storage_manager_.open_array(result_array_schema);
    
  // Get the tile ids that overlap with the range
  std::vector<std::pair<uint64_t, bool> > overlapping_tile_ids;
  storage_manager_.get_overlapping_tile_ids(ad, range, &overlapping_tile_ids);

  // Initialize tile iterators
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it = 
      overlapping_tile_ids.begin();
  std::vector<std::pair<uint64_t, bool> >::const_iterator tile_id_it_end =
      overlapping_tile_ids.end();


  // Auxiliary variable storing the number of skipped cells when investigating
  // a tile partially overlapping the range. It is used to advance only the
  // coordinates iterator when a cell is not in range, and then efficiently
  // advance the attribute iterators only when a cell falls in the range.
  int64_t skipped;

  // Iterate over all overlapping tiles
  for(; tile_id_it != tile_id_it_end; ++tile_id_it) {
    // Create result tiles and load input array tiles 
    new_tiles(result_array_schema, tile_id_it->first, result_tiles); 
    get_tiles(ad, tile_id_it->first, tiles); 
    initialize_cell_its(tiles, attribute_num, cell_its, cell_it_end); 
    skipped = 0;
 
    if(tile_id_it->second) { // Full overlap
      while(cell_its[attribute_num] != cell_it_end) {
        append_cell(cell_its, result_tiles, attribute_num);
        advance_cell_its(attribute_num, cell_its);
      }
    } else { // Partial overlap
      while(cell_its[attribute_num] != cell_it_end) {
        if(cell_its[attribute_num].cell_inside_range(range)) {
          advance_cell_its(attribute_num, cell_its, skipped);
          skipped = 0;
          append_cell(cell_its, result_tiles, attribute_num);
          advance_cell_its(attribute_num, cell_its);
        } else { // Advance only the coordinates cell iterator
          ++skipped;
          ++cell_its[attribute_num];
        }
      }
    }
      
    // Send new tiles to storage manager
    store_tiles(result_ad, result_tiles);
  } 
  
  // Close result array
  storage_manager_.close_array(result_ad);

  // Clean up 
  delete [] tiles;
  delete [] result_tiles;
  delete [] cell_its;
}

