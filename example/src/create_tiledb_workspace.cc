#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include "c_api.h"

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Needs 1 argument <workspace_directory>\n";
    exit(-1);
  }
  auto workspace = argv[1];
  //Create workspace if it does not exist
  struct stat st;
  auto status = stat(workspace, &st);
  int returnval = 0;
  //Exists and is not a directory
  if(status >= 0 && !S_ISDIR(st.st_mode))
  {
    std::cerr << "Workspace path " << workspace << " exists and is not a directory\n";
    returnval = -1;
  }
  else
  {
    if(status >= 0)
      std::cerr << "Directory " << workspace << " exists - doing nothing\n";
    else  //Doesn't exist, create workspace
    {
      TileDB_CTX* tiledb_ctx = 0;
      /*Initialize context with default params*/
      tiledb_ctx_init(&tiledb_ctx, NULL);
      if(tiledb_workspace_create(tiledb_ctx, workspace) != TILEDB_OK)
      {
        std::cerr << "Failed to create workspace "<<workspace<<"\n";
        returnval = -1;
      }
      else
        std::cerr << "Created workspace "<<workspace<<"\n";
      tiledb_ctx_finalize(tiledb_ctx);
      free(tiledb_ctx);
    }
  }
  return returnval;
}
