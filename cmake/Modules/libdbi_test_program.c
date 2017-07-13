#include <stdio.h>
#include <dbi/dbi.h>

int main() {
  dbi_conn conn;
  dbi_result result;
  dbi_inst instance;

  dbi_initialize_r(NULL, &instance);
  conn = dbi_conn_new_r("pgsql", instance);

  dbi_conn_set_option(conn, "host", "localhost");

  if (dbi_conn_connect(conn) < 0)
    exit(0);
  else
  {
    result = dbi_conn_query(conn, "SELECT * from table");
    if (result) {
      while (dbi_result_next_row(result)) {
        unsigned idnumber = dbi_result_get_uint(result, "id");
      }
      dbi_result_free(result);
    }
    dbi_conn_close(conn);
  }
  dbi_shutdown_r(instance);
  return 0;
}
