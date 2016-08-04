namespace java tutorial

typedef i32 int
typedef i64 int64

service GenomicsDBService {
    int getVariantContext(1:String loaderJson, 2:String queryJson),
}