--
-- PostgreSQL database dump
--

-- Dumped from database version 9.5.7
-- Dumped by pg_dump version 9.5.7

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET search_path = public, pg_catalog;

--
-- Name: field_combine_optype; Type: TYPE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TYPE field_combine_optype AS ENUM (
    'sum',
    'mean',
    'median',
    'move_to_FORMAT',
    'element_wise_sum',
    'concatenate'
);


ALTER TYPE field_combine_optype OWNER TO karthikg_genomicsdb;

--
-- Name: length_enum; Type: TYPE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TYPE length_enum AS ENUM (
    'A',
    'R',
    'G',
    'VAR',
    'NUM'
);


ALTER TYPE length_enum OWNER TO karthikg_genomicsdb;

--
-- Name: type_enum; Type: TYPE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TYPE type_enum AS ENUM (
    'Integer',
    'String',
    'Float',
    'Flag'
);


ALTER TYPE type_enum OWNER TO karthikg_genomicsdb;

--
-- Name: increment_next_column_in_reference_set_pgsql(); Type: FUNCTION; Schema: public; Owner: karthikg_genomicsdb
--

CREATE FUNCTION increment_next_column_in_reference_set_pgsql() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
    DECLARE
        updated_next_tiledb_column_offset bigint;
        padded_reference_length bigint;
    BEGIN
        padded_reference_length = CAST( CAST(NEW.length AS DOUBLE PRECISION)*(SELECT tiledb_reference_offset_padding_factor FROM reference_set WHERE id=NEW.reference_set_id) AS BIGINT);
        UPDATE reference_set SET next_tiledb_column_offset=
            CASE
                WHEN NEW.tiledb_column_offset IS NULL THEN next_tiledb_column_offset+padded_reference_length
                WHEN NEW.tiledb_column_offset+padded_reference_length>next_tiledb_column_offset THEN NEW.tiledb_column_offset+padded_reference_length
                ELSE next_tiledb_column_offset
            END
        WHERE id = NEW.reference_set_id RETURNING next_tiledb_column_offset INTO updated_next_tiledb_column_offset;
        IF NEW.tiledb_column_offset IS NULL THEN
            NEW.tiledb_column_offset = updated_next_tiledb_column_offset-padded_reference_length;
        END IF;
        RETURN NEW;
    END;
    $$;


ALTER FUNCTION public.increment_next_column_in_reference_set_pgsql() OWNER TO karthikg_genomicsdb;

--
-- Name: increment_num_rows_in_db_array_pgsql(); Type: FUNCTION; Schema: public; Owner: karthikg_genomicsdb
--

CREATE FUNCTION increment_num_rows_in_db_array_pgsql() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
    DECLARE
        updated_num_rows bigint;
    BEGIN
        UPDATE db_array SET num_rows=
            CASE
               WHEN NEW.tile_row_id IS NULL THEN num_rows+1
               WHEN NEW.tile_row_id >= num_rows THEN NEW.tile_row_id+1
               ELSE num_rows
            END
        WHERE id=NEW.db_array_id RETURNING num_rows INTO updated_num_rows;
        IF NEW.tile_row_id IS NULL THEN
            NEW.tile_row_id = updated_num_rows-1;
        END IF;
        RETURN NEW;
    END;
    $$;


ALTER FUNCTION public.increment_num_rows_in_db_array_pgsql() OWNER TO karthikg_genomicsdb;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: alembic_version; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE alembic_version (
    version_num character varying(32) NOT NULL
);


ALTER TABLE alembic_version OWNER TO karthikg_genomicsdb;

--
-- Name: callset_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE callset_id_seq
    START WITH 0
    INCREMENT BY 1
    MINVALUE 0
    NO MAXVALUE
    CACHE 1;


ALTER TABLE callset_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: callset; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE callset (
    id bigint DEFAULT nextval('callset_id_seq'::regclass) NOT NULL,
    guid character varying(36) NOT NULL,
    name text NOT NULL,
    created bigint NOT NULL,
    updated bigint NOT NULL,
    info bytea,
    source_sample_id bigint NOT NULL,
    target_sample_id bigint NOT NULL
);


ALTER TABLE callset OWNER TO karthikg_genomicsdb;

--
-- Name: callset_to_db_array_association; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE callset_to_db_array_association (
    callset_id bigint NOT NULL,
    db_array_id bigint NOT NULL,
    tile_row_id bigint
);


ALTER TABLE callset_to_db_array_association OWNER TO karthikg_genomicsdb;

--
-- Name: callset_variant_set; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE callset_variant_set (
    callset_id bigint NOT NULL,
    variant_set_id bigint NOT NULL
);


ALTER TABLE callset_variant_set OWNER TO karthikg_genomicsdb;

--
-- Name: db_array; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE db_array (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    reference_set_id bigint NOT NULL,
    workspace_id bigint NOT NULL,
    name text NOT NULL,
    num_rows bigint DEFAULT 0 NOT NULL,
    field_set_id bigint NOT NULL
);


ALTER TABLE db_array OWNER TO karthikg_genomicsdb;

--
-- Name: db_array_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE db_array_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE db_array_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: db_array_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE db_array_id_seq OWNED BY db_array.id;


--
-- Name: db_row_tile_row_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE db_row_tile_row_id_seq
    START WITH 0
    INCREMENT BY 1
    MINVALUE 0
    NO MAXVALUE
    CACHE 1;


ALTER TABLE db_row_tile_row_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: field; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE field (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    field character varying(32) NOT NULL,
    field_set_id bigint NOT NULL,
    type type_enum,
    is_filter boolean NOT NULL,
    is_format boolean NOT NULL,
    is_info boolean NOT NULL,
    length_type length_enum,
    length_intval integer DEFAULT 1,
    field_combine_op field_combine_optype
);


ALTER TABLE field OWNER TO karthikg_genomicsdb;

--
-- Name: field_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE field_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE field_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: field_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE field_id_seq OWNED BY field.id;


--
-- Name: field_set; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE field_set (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    description text
);


ALTER TABLE field_set OWNER TO karthikg_genomicsdb;

--
-- Name: field_set_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE field_set_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE field_set_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: field_set_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE field_set_id_seq OWNED BY field_set.id;


--
-- Name: individual; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE individual (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    name text NOT NULL,
    info bytea,
    record_create_time text,
    record_update_time text
);


ALTER TABLE individual OWNER TO karthikg_genomicsdb;

--
-- Name: individual_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE individual_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE individual_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: individual_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE individual_id_seq OWNED BY individual.id;


--
-- Name: reference; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE reference (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    length bigint NOT NULL,
    reference_set_id bigint NOT NULL,
    md5_checksum character varying(32),
    name text NOT NULL,
    source_uri text,
    is_derived boolean,
    source_divergence double precision,
    ncbi_taxon_id integer,
    tiledb_column_offset bigint
);


ALTER TABLE reference OWNER TO karthikg_genomicsdb;

--
-- Name: reference_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE reference_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE reference_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: reference_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE reference_id_seq OWNED BY reference.id;


--
-- Name: reference_set; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE reference_set (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    md5_checksum character varying(32),
    description text,
    source_uri text,
    is_derived boolean,
    ncbi_taxon_id integer,
    assembly_id character varying(100),
    next_tiledb_column_offset bigint DEFAULT 0 NOT NULL,
    tiledb_reference_offset_padding_factor double precision DEFAULT 1.10 NOT NULL
);


ALTER TABLE reference_set OWNER TO karthikg_genomicsdb;

--
-- Name: reference_set_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE reference_set_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE reference_set_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: reference_set_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE reference_set_id_seq OWNED BY reference_set.id;


--
-- Name: reference_set_source_accession; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE reference_set_source_accession (
    reference_set_id bigint NOT NULL,
    source_accession_id bigint NOT NULL
);


ALTER TABLE reference_set_source_accession OWNER TO karthikg_genomicsdb;

--
-- Name: reference_source_accession; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE reference_source_accession (
    reference_id bigint NOT NULL,
    source_accession_id bigint NOT NULL
);


ALTER TABLE reference_source_accession OWNER TO karthikg_genomicsdb;

--
-- Name: sample; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE sample (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    individual_id bigint NOT NULL,
    name text NOT NULL,
    info text
);


ALTER TABLE sample OWNER TO karthikg_genomicsdb;

--
-- Name: sample_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE sample_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE sample_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: sample_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE sample_id_seq OWNED BY sample.id;


--
-- Name: source_accession; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE source_accession (
    id bigint NOT NULL,
    accession_id text NOT NULL
);


ALTER TABLE source_accession OWNER TO karthikg_genomicsdb;

--
-- Name: source_accession_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE source_accession_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE source_accession_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: source_accession_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE source_accession_id_seq OWNED BY source_accession.id;


--
-- Name: variant_set; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE variant_set (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    name text,
    reference_set_id bigint NOT NULL,
    dataset_id text,
    variant_set_metadata text
);


ALTER TABLE variant_set OWNER TO karthikg_genomicsdb;

--
-- Name: variant_set_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE variant_set_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE variant_set_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: variant_set_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE variant_set_id_seq OWNED BY variant_set.id;


--
-- Name: workspace; Type: TABLE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TABLE workspace (
    id bigint NOT NULL,
    guid character varying(36) NOT NULL,
    name text NOT NULL
);


ALTER TABLE workspace OWNER TO karthikg_genomicsdb;

--
-- Name: workspace_id_seq; Type: SEQUENCE; Schema: public; Owner: karthikg_genomicsdb
--

CREATE SEQUENCE workspace_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE workspace_id_seq OWNER TO karthikg_genomicsdb;

--
-- Name: workspace_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: karthikg_genomicsdb
--

ALTER SEQUENCE workspace_id_seq OWNED BY workspace.id;


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY db_array ALTER COLUMN id SET DEFAULT nextval('db_array_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field ALTER COLUMN id SET DEFAULT nextval('field_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field_set ALTER COLUMN id SET DEFAULT nextval('field_set_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY individual ALTER COLUMN id SET DEFAULT nextval('individual_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference ALTER COLUMN id SET DEFAULT nextval('reference_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_set ALTER COLUMN id SET DEFAULT nextval('reference_set_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY sample ALTER COLUMN id SET DEFAULT nextval('sample_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY source_accession ALTER COLUMN id SET DEFAULT nextval('source_accession_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY variant_set ALTER COLUMN id SET DEFAULT nextval('variant_set_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY workspace ALTER COLUMN id SET DEFAULT nextval('workspace_id_seq'::regclass);


--
-- Data for Name: alembic_version; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY alembic_version (version_num) FROM stdin;
4f93dc7aa8e8
\.


--
-- Data for Name: callset; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY callset (id, guid, name, created, updated, info, source_sample_id, target_sample_id) FROM stdin;
0	d6bc50b2-0d2a-499b-a778-74b56c328e9d	HG01958	1499735659127	1499735659127	\N	1	1
1	cb2f7d4e-7c3c-49e8-b7e0-27b34cf8cdb3	HG01530	1499735659147	1499735659147	\N	2	2
2	522374db-af06-4907-a799-e518f9bae27f	HG00141	1499735659151	1499735659151	\N	3	3
\.


--
-- Name: callset_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('callset_id_seq', 2, true);


--
-- Data for Name: callset_to_db_array_association; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY callset_to_db_array_association (callset_id, db_array_id, tile_row_id) FROM stdin;
2	1	0
0	1	1
1	1	2
\.


--
-- Data for Name: callset_variant_set; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY callset_variant_set (callset_id, variant_set_id) FROM stdin;
0	1
1	1
2	1
\.


--
-- Data for Name: db_array; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY db_array (id, guid, reference_set_id, workspace_id, name, num_rows, field_set_id) FROM stdin;
1	cafbfca8-4ce6-4b33-bdef-eaceca700f09	1	1	t0_1_2	3	1
\.


--
-- Name: db_array_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('db_array_id_seq', 1, true);


--
-- Name: db_row_tile_row_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('db_row_tile_row_id_seq', 0, false);


--
-- Data for Name: field; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY field (id, guid, field, field_set_id, type, is_filter, is_format, is_info, length_type, length_intval, field_combine_op) FROM stdin;
1	deb419bd-79b8-4c4f-9e4e-4ff1de8e1dcc	BaseQRankSum	1	Float	f	f	t	NUM	1	\N
2	7f1ededb-5d24-4aaf-bbd9-f84d3143e4e1	MQRankSum	1	Float	f	f	t	NUM	1	\N
4	bebc66ae-89f7-4bd0-bb19-ba074ad43c01	DP	1	Integer	f	t	t	NUM	1	\N
5	a5e5b015-6fb3-4c38-b694-fc19976dcda0	MLEAC	1	Integer	f	f	t	A	1	\N
6	b6a6a120-1263-4290-a2c1-f0b3521fac2a	MLEAF	1	Float	f	f	t	A	1	\N
9	f17b4b82-0003-4d17-ba6b-a1822843f632	HaplotypeScore	1	Float	f	f	t	NUM	1	\N
10	7517618e-8899-4b1a-850f-0e645c9927c8	PL	1	Integer	f	t	f	G	1	\N
12	ea3122af-cf50-43f4-89ad-c9558cb0d13c	END	1	Integer	f	f	t	NUM	1	\N
14	57193118-27e5-468f-a183-14543bea4907	GQ	1	Integer	f	t	f	NUM	1	\N
15	c271a1e9-18a4-45d7-b18a-26a59673ff70	MQ0	1	Integer	f	f	t	NUM	1	\N
16	6b6ff2da-c75b-4902-9205-710576fd68ce	RAW_MQ	1	Float	f	f	t	NUM	1	\N
17	1ff8623e-c9d3-4fe3-be70-03c30dddb718	PASS	1	\N	t	f	f	\N	1	\N
18	7c8f0c35-c379-40bb-a3cc-a3269bfbf29e	ReadPosRankSum	1	Float	f	f	t	NUM	1	\N
19	7699415b-4b87-402a-bbdb-1db2e03f971d	MIN_DP	1	Integer	f	t	f	NUM	1	\N
20	d33dc19c-6516-422d-aceb-4509669f43ab	LowQual	1	\N	t	f	f	\N	1	\N
21	528e0d7d-b868-4d39-be0a-d834701706f7	InbreedingCoeff	1	Float	f	f	t	NUM	1	\N
22	18a8ded9-dddc-45cf-be8a-5972211f2f3a	MQ	1	Float	f	f	t	NUM	1	\N
23	aab4c86a-6745-4a76-842d-6ad57f127e12	ClippingRankSum	1	Float	f	f	t	NUM	1	\N
24	50b303f9-84d8-4871-9bb4-159b830da182	SB	1	Integer	f	t	f	NUM	4	\N
7	d2e99501-d80c-4cb7-9fd2-002539b6a9f6	PID	1	String	f	t	f	VAR	1	\N
8	d43651ff-cc74-4beb-9297-9889c807bee1	PGT	1	String	f	t	f	VAR	1	\N
11	9397f903-fb5e-4f4f-94b3-097f2b804fa7	GT	1	Integer	f	t	f	VAR	1	\N
13	c373f452-b1f4-4812-8fc2-820cbd28ffea	AD	1	Integer	f	t	f	R	1	\N
3	889e81a8-cda2-42fc-b6ba-6b6032b1c0fa	DS	1	Flag	f	f	t	NUM	1	\N
25	daae5fad-c32d-4fc6-800e-3cc34e2376b2	ID	1	String	f	f	f	VAR	1	\N
\.


--
-- Name: field_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('field_id_seq', 25, true);


--
-- Data for Name: field_set; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY field_set (id, guid, description) FROM stdin;
1	8b058c58-bb52-4f1c-b440-21c722553f01	hg19
\.


--
-- Name: field_set_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('field_set_id_seq', 1, true);


--
-- Data for Name: individual; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY individual (id, guid, name, info, record_create_time, record_update_time) FROM stdin;
1	61011f96-afec-4b0c-9f12-dd3358434e30	Individual_HG01958	\N	2017-07-10 18:14:1919.191919	2017-07-10 18:14:1919.191919
2	e6f59f65-b5e3-4969-b24b-77d3b78dc653	Individual_HG01530	\N	2017-07-10 18:14:1919.191919	2017-07-10 18:14:1919.191919
3	9fe5dad6-9b12-48c9-b0dd-98b002fa22d1	Individual_HG00141	\N	2017-07-10 18:14:1919.191919	2017-07-10 18:14:1919.191919
\.


--
-- Name: individual_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('individual_id_seq', 3, true);


--
-- Data for Name: reference; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY reference (id, guid, length, reference_set_id, md5_checksum, name, source_uri, is_derived, source_divergence, ncbi_taxon_id, tiledb_column_offset) FROM stdin;
1	7dc05155-d771-4d30-ae6a-6a6ec508d036	249250621	1	\N	1	\N	\N	\N	\N	0
2	285024a1-b206-4922-a70e-8c0aaadba006	243199373	1	\N	2	\N	\N	\N	\N	274175683
3	2f7faefb-bc4b-48ea-8d02-86f89e4fddf6	198022430	1	\N	3	\N	\N	\N	\N	541694993
4	7bf814f9-cc8a-4790-a150-80c4181d827f	191154276	1	\N	4	\N	\N	\N	\N	759519666
5	0bb7d6d0-9e40-4379-bab0-ff28b06b4977	180915260	1	\N	5	\N	\N	\N	\N	969789370
6	f0dea69c-a90d-468e-bfb2-43241415eb3a	171115067	1	\N	6	\N	\N	\N	\N	1168796156
7	5ccc3ed9-214d-4b82-9af4-5a568a778da1	159138663	1	\N	7	\N	\N	\N	\N	1357022730
8	b70c1738-a742-45d9-9288-094f44a0c064	146364022	1	\N	8	\N	\N	\N	\N	1532075259
9	e94850f1-a575-4d95-bd22-b91ba073549c	141213431	1	\N	9	\N	\N	\N	\N	1693075683
10	f570f3e5-8163-4d7b-9f6f-84748fc42f16	135534747	1	\N	10	\N	\N	\N	\N	1848410457
11	acd49448-8daf-4c15-b6dc-d6f15fbe6897	135006516	1	\N	11	\N	\N	\N	\N	1997498679
12	e47ebdb8-f77a-4100-94c7-120d604eb21a	133851895	1	\N	12	\N	\N	\N	\N	2146005847
13	b98e81d2-bb1b-4eb3-a33d-fba058b29d83	115169878	1	\N	13	\N	\N	\N	\N	2293242931
14	48dffafb-5f5e-4238-a2a7-83479d4b79a1	107349540	1	\N	14	\N	\N	\N	\N	2419929797
15	7997855e-69f1-4a0c-af79-224f7f086d27	102531392	1	\N	15	\N	\N	\N	\N	2538014291
16	c22e6a78-4c47-43fb-954b-6734065fcd6a	90354753	1	\N	16	\N	\N	\N	\N	2650798822
17	18d3fbef-338a-49fb-945c-06eb381fa851	81195210	1	\N	17	\N	\N	\N	\N	2750189050
18	afcde2f7-777a-402e-9c1a-5cae5daff82b	78077248	1	\N	18	\N	\N	\N	\N	2839503781
19	c633fb4e-96bd-44ec-8fef-46e060411226	59128983	1	\N	19	\N	\N	\N	\N	2925388754
20	13e4b3fd-11c4-4cc2-91f1-ac1c4cacab85	63025520	1	\N	20	\N	\N	\N	\N	2990430635
21	699664f9-d4ec-48ff-a8f6-fcf188779265	48129895	1	\N	21	\N	\N	\N	\N	3059758707
22	db0e851e-ab48-4817-935b-acc93657051f	51304566	1	\N	22	\N	\N	\N	\N	3112701592
23	1b1609e3-7fe0-48a9-a226-19503d7fdc71	155270560	1	\N	X	\N	\N	\N	\N	3169136615
24	a3285396-79c8-43a2-9f94-3f1e3ef5c38c	59373566	1	\N	Y	\N	\N	\N	\N	3339934231
25	7c5c624c-99ca-4ad4-aaf8-e39575b66195	4262	1	\N	GL000207.1	\N	\N	\N	\N	3405245154
26	c044a1c9-4148-473e-a515-ce287ec2f729	15008	1	\N	GL000226.1	\N	\N	\N	\N	3405249842
27	9f8e2924-69e4-4a12-99b1-207d1cd1f550	19913	1	\N	GL000229.1	\N	\N	\N	\N	3405266351
28	63e68d31-7d41-4f18-9729-e17496fdf362	27386	1	\N	GL000231.1	\N	\N	\N	\N	3405288255
29	cdc28cfc-437d-4c8b-9358-31147a5c4ad6	27682	1	\N	GL000210.1	\N	\N	\N	\N	3405318380
30	d32f411f-97a9-4f25-9fbc-d5ad06bd3dd3	33824	1	\N	GL000239.1	\N	\N	\N	\N	3405348830
31	30be12ef-7eaa-40b5-98ca-c18ee20c02ed	34474	1	\N	GL000235.1	\N	\N	\N	\N	3405386036
32	c798c6d8-0314-46af-9de7-a609e2aecfec	36148	1	\N	GL000201.1	\N	\N	\N	\N	3405423957
33	308668a0-e790-4784-af6b-9b6bd5335d83	36422	1	\N	GL000247.1	\N	\N	\N	\N	3405463720
34	63bae760-fd1a-48b1-aafc-d0391c691e33	36651	1	\N	GL000245.1	\N	\N	\N	\N	3405503784
35	7b57fa14-df36-43ac-b918-6428ada96c27	37175	1	\N	GL000197.1	\N	\N	\N	\N	3405544100
36	6b4d068b-9e0c-4b47-93f7-1e0f35ffd847	37498	1	\N	GL000203.1	\N	\N	\N	\N	3405584992
37	072807ea-f7a7-453b-b03d-1b7d049cc4a8	38154	1	\N	GL000246.1	\N	\N	\N	\N	3405626240
38	28e515a5-e13d-4e7b-bc95-951737a95e41	38502	1	\N	GL000249.1	\N	\N	\N	\N	3405668209
39	ab1e36c2-5483-4615-bbcc-db693d14947f	38914	1	\N	GL000196.1	\N	\N	\N	\N	3405710561
40	6882f73a-9b3d-4a86-91f6-e491825c9d27	39786	1	\N	GL000248.1	\N	\N	\N	\N	3405753366
41	8cd1d37c-99c5-4e36-a26e-47835f36f147	39929	1	\N	GL000244.1	\N	\N	\N	\N	3405797131
42	a0fb4050-92f5-48a8-af24-0f66a220f51b	39939	1	\N	GL000238.1	\N	\N	\N	\N	3405841053
43	40473b0e-6886-479f-b6ed-7c14ae099222	40103	1	\N	GL000202.1	\N	\N	\N	\N	3405884986
44	de410c14-70ab-44f1-bb59-2487059e6532	40531	1	\N	GL000234.1	\N	\N	\N	\N	3405929099
45	7a9ca471-debb-4894-b63c-e0fb3ad6b081	40652	1	\N	GL000232.1	\N	\N	\N	\N	3405973683
46	407f6d24-f93d-4bf5-8655-bba98db2ed88	41001	1	\N	GL000206.1	\N	\N	\N	\N	3406018400
47	344f1957-c1a6-4818-bfb2-e32bb8b4a204	41933	1	\N	GL000240.1	\N	\N	\N	\N	3406063501
48	df1532c5-72bf-4fa0-bec9-281eac19036b	41934	1	\N	GL000236.1	\N	\N	\N	\N	3406109627
49	231432e2-38a9-4b31-9f1f-ef4fea23e432	42152	1	\N	GL000241.1	\N	\N	\N	\N	3406155754
50	20bc97d7-a20f-429e-8a52-89ebcbd4563c	43341	1	\N	GL000243.1	\N	\N	\N	\N	3406202121
51	93436176-1643-49dd-afde-58f9326bd93e	43523	1	\N	GL000242.1	\N	\N	\N	\N	3406249796
52	2973e143-96e4-4d2d-ad99-f599312c99e9	43691	1	\N	GL000230.1	\N	\N	\N	\N	3406297671
53	c23dad1b-010c-4631-a34f-b7b010fe359c	45867	1	\N	GL000237.1	\N	\N	\N	\N	3406345731
54	408bd8cd-eb84-426a-b7ac-52a2af965e9d	45941	1	\N	GL000233.1	\N	\N	\N	\N	3406396185
55	6754285c-ae1d-4c3e-9f8a-9b43ceb4f1b1	81310	1	\N	GL000204.1	\N	\N	\N	\N	3406446720
56	008c416e-532a-4100-8ce3-ee49994c25da	90085	1	\N	GL000198.1	\N	\N	\N	\N	3406536161
57	98e58d9b-d4e9-4c44-a057-8903100e0f32	92689	1	\N	GL000208.1	\N	\N	\N	\N	3406635255
58	6298c60f-26c8-4ef7-aa68-afe0a5f19308	106433	1	\N	GL000191.1	\N	\N	\N	\N	3406737213
59	1c28d69b-563b-4520-93f8-35ae66c5d98c	128374	1	\N	GL000227.1	\N	\N	\N	\N	3406854289
60	1a9b0212-0e38-47f7-b492-c8fbc5582da8	129120	1	\N	GL000228.1	\N	\N	\N	\N	3406995500
61	6ceb9863-d37a-417c-9714-8a43c8dda342	137718	1	\N	GL000214.1	\N	\N	\N	\N	3407137532
62	10eaa66d-4808-45eb-9544-bc772f6d73ff	155397	1	\N	GL000221.1	\N	\N	\N	\N	3407289022
63	3e0820d8-e7c9-4a66-afc7-06536bd69806	159169	1	\N	GL000209.1	\N	\N	\N	\N	3407459959
64	9015acd4-64f9-49e3-80dc-150d36902279	161147	1	\N	GL000218.1	\N	\N	\N	\N	3407635045
65	5688645a-2d1f-48bd-a59e-5dee875dcc0b	161802	1	\N	GL000220.1	\N	\N	\N	\N	3407812307
66	e841db9f-4234-4248-9e5c-314385570264	164239	1	\N	GL000213.1	\N	\N	\N	\N	3407990289
67	82b9dae4-a87d-4136-93f0-9fa9470a9f9d	166566	1	\N	GL000211.1	\N	\N	\N	\N	3408170952
68	a29ab00d-1462-4849-a742-814615de809a	169874	1	\N	GL000199.1	\N	\N	\N	\N	3408354175
69	64c87c20-b6d7-4083-99af-cdd5330a41f1	172149	1	\N	GL000217.1	\N	\N	\N	\N	3408541036
70	0e7b53f2-3099-4067-8449-1c85775f5b48	172294	1	\N	GL000216.1	\N	\N	\N	\N	3408730400
71	fd7e1caa-fddb-4e0a-a521-0df90f65a619	172545	1	\N	GL000215.1	\N	\N	\N	\N	3408919923
72	2911d4a8-d3b3-49fd-82ca-aefbbcadc8a5	174588	1	\N	GL000205.1	\N	\N	\N	\N	3409109723
73	4983b58c-603d-4758-925a-0267da26101c	179198	1	\N	GL000219.1	\N	\N	\N	\N	3409301770
74	d2529bdd-a5a2-4e45-ab1b-c7e9a5553d54	179693	1	\N	GL000224.1	\N	\N	\N	\N	3409498888
75	ff9fab1c-8baf-4385-a813-748f76ebf540	180455	1	\N	GL000223.1	\N	\N	\N	\N	3409696550
76	143a67a9-de0f-4986-991f-cb08ee4db0c6	182896	1	\N	GL000195.1	\N	\N	\N	\N	3409895051
77	967529cf-ac1a-49cd-a4cd-686275be47f1	186858	1	\N	GL000212.1	\N	\N	\N	\N	3410096237
78	626dbe1a-14ea-4d66-a15f-3fbc564b47fd	186861	1	\N	GL000222.1	\N	\N	\N	\N	3410301781
79	17cf5b6c-0067-415f-9724-7254476dcc4a	187035	1	\N	GL000200.1	\N	\N	\N	\N	3410507328
80	984c7955-031e-499a-ba63-92c9ea9377df	189789	1	\N	GL000193.1	\N	\N	\N	\N	3410713067
81	9acb6c01-7545-4274-b4d1-c5d33a9ad3f5	191469	1	\N	GL000194.1	\N	\N	\N	\N	3410921835
82	6ae207c6-7a55-4109-b25d-3e6e5bc8cb84	211173	1	\N	GL000225.1	\N	\N	\N	\N	3411132451
83	858be628-7669-4b03-9676-03ca94b7b01c	547496	1	\N	GL000192.1	\N	\N	\N	\N	3411364741
84	97fd2c29-d7f0-40c7-8f38-5993dd6cdf46	171823	1	\N	NC_007605	\N	\N	\N	\N	3411966987
85	c13d3d9f-e00f-4b85-acf1-fd41b123d336	16569	1	\N	M	\N	\N	\N	\N	3412155992
\.


--
-- Name: reference_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('reference_id_seq', 85, true);


--
-- Data for Name: reference_set; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY reference_set (id, guid, md5_checksum, description, source_uri, is_derived, ncbi_taxon_id, assembly_id, next_tiledb_column_offset, tiledb_reference_offset_padding_factor) FROM stdin;
1	0831f68d-23f6-46d6-a79f-f6abd1729f8c	\N	\N	\N	\N	\N	hg19	3412174218	1.10000000000000009
\.


--
-- Name: reference_set_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('reference_set_id_seq', 1, true);


--
-- Data for Name: reference_set_source_accession; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY reference_set_source_accession (reference_set_id, source_accession_id) FROM stdin;
\.


--
-- Data for Name: reference_source_accession; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY reference_source_accession (reference_id, source_accession_id) FROM stdin;
\.


--
-- Data for Name: sample; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY sample (id, guid, individual_id, name, info) FROM stdin;
1	a8cf8937-825a-4589-8729-f0e7f17d92f0	1	HG01958	{"type": "source"}
2	ea62aa59-2318-44ec-a212-5167abfe52fd	2	HG01530	{"type": "source"}
3	15d0fa74-455c-4922-b62f-77398fc70c4c	3	HG00141	{"type": "source"}
\.


--
-- Name: sample_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('sample_id_seq', 3, true);


--
-- Data for Name: source_accession; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY source_accession (id, accession_id) FROM stdin;
\.


--
-- Name: source_accession_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('source_accession_id_seq', 1, false);


--
-- Data for Name: variant_set; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY variant_set (id, guid, name, reference_set_id, dataset_id, variant_set_metadata) FROM stdin;
1	f1b99e17-543b-4797-96a6-ce85c15082e7	\N	1	ws	\N
\.


--
-- Name: variant_set_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('variant_set_id_seq', 1, true);


--
-- Data for Name: workspace; Type: TABLE DATA; Schema: public; Owner: karthikg_genomicsdb
--

COPY workspace (id, guid, name) FROM stdin;
1	44400e0a-5095-4cfa-a577-b44b4b58d071	/tmp/ws
\.


--
-- Name: workspace_id_seq; Type: SEQUENCE SET; Schema: public; Owner: karthikg_genomicsdb
--

SELECT pg_catalog.setval('workspace_id_seq', 1, true);


--
-- Name: callset_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset
    ADD CONSTRAINT callset_guid_key UNIQUE (guid);


--
-- Name: callset_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset
    ADD CONSTRAINT callset_pkey PRIMARY KEY (id);


--
-- Name: callset_variant_set_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset_variant_set
    ADD CONSTRAINT callset_variant_set_pkey PRIMARY KEY (callset_id, variant_set_id);


--
-- Name: db_array_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY db_array
    ADD CONSTRAINT db_array_guid_key UNIQUE (guid);


--
-- Name: db_array_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY db_array
    ADD CONSTRAINT db_array_pkey PRIMARY KEY (id);


--
-- Name: field_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field
    ADD CONSTRAINT field_guid_key UNIQUE (guid);


--
-- Name: field_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field
    ADD CONSTRAINT field_pkey PRIMARY KEY (id);


--
-- Name: field_set_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field_set
    ADD CONSTRAINT field_set_guid_key UNIQUE (guid);


--
-- Name: field_set_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field_set
    ADD CONSTRAINT field_set_pkey PRIMARY KEY (id);


--
-- Name: individual_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY individual
    ADD CONSTRAINT individual_guid_key UNIQUE (guid);


--
-- Name: individual_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY individual
    ADD CONSTRAINT individual_pkey PRIMARY KEY (id);


--
-- Name: primary_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset_to_db_array_association
    ADD CONSTRAINT primary_key PRIMARY KEY (callset_id, db_array_id);


--
-- Name: reference_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference
    ADD CONSTRAINT reference_guid_key UNIQUE (guid);


--
-- Name: reference_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference
    ADD CONSTRAINT reference_pkey PRIMARY KEY (id);


--
-- Name: reference_set_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_set
    ADD CONSTRAINT reference_set_guid_key UNIQUE (guid);


--
-- Name: reference_set_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_set
    ADD CONSTRAINT reference_set_pkey PRIMARY KEY (id);


--
-- Name: reference_set_source_accession_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_set_source_accession
    ADD CONSTRAINT reference_set_source_accession_pkey PRIMARY KEY (reference_set_id, source_accession_id);


--
-- Name: reference_source_accession_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_source_accession
    ADD CONSTRAINT reference_source_accession_pkey PRIMARY KEY (reference_id, source_accession_id);


--
-- Name: sample_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_guid_key UNIQUE (guid);


--
-- Name: sample_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_pkey PRIMARY KEY (id);


--
-- Name: source_accession_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY source_accession
    ADD CONSTRAINT source_accession_pkey PRIMARY KEY (id);


--
-- Name: unique_name_per_field_set_constraint; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field
    ADD CONSTRAINT unique_name_per_field_set_constraint UNIQUE (field_set_id, field);


--
-- Name: unique_name_per_reference_set_constraint; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference
    ADD CONSTRAINT unique_name_per_reference_set_constraint UNIQUE (reference_set_id, name);


--
-- Name: variant_set_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY variant_set
    ADD CONSTRAINT variant_set_guid_key UNIQUE (guid);


--
-- Name: variant_set_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY variant_set
    ADD CONSTRAINT variant_set_pkey PRIMARY KEY (id);


--
-- Name: workspace_guid_key; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY workspace
    ADD CONSTRAINT workspace_guid_key UNIQUE (guid);


--
-- Name: workspace_pkey; Type: CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY workspace
    ADD CONSTRAINT workspace_pkey PRIMARY KEY (id);


--
-- Name: db_array_id_tile_row_id_idx; Type: INDEX; Schema: public; Owner: karthikg_genomicsdb
--

CREATE UNIQUE INDEX db_array_id_tile_row_id_idx ON callset_to_db_array_association USING btree (db_array_id, tile_row_id);


--
-- Name: unique_reference_set_id_offset_idx; Type: INDEX; Schema: public; Owner: karthikg_genomicsdb
--

CREATE UNIQUE INDEX unique_reference_set_id_offset_idx ON reference USING btree (reference_set_id, tiledb_column_offset);


--
-- Name: increment_next_column_in_reference_set; Type: TRIGGER; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TRIGGER increment_next_column_in_reference_set BEFORE INSERT ON reference FOR EACH ROW EXECUTE PROCEDURE increment_next_column_in_reference_set_pgsql();


--
-- Name: increment_num_rows_in_db_array; Type: TRIGGER; Schema: public; Owner: karthikg_genomicsdb
--

CREATE TRIGGER increment_num_rows_in_db_array BEFORE INSERT ON callset_to_db_array_association FOR EACH ROW EXECUTE PROCEDURE increment_num_rows_in_db_array_pgsql();


--
-- Name: callset_source_sample_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset
    ADD CONSTRAINT callset_source_sample_id_fkey FOREIGN KEY (source_sample_id) REFERENCES sample(id);


--
-- Name: callset_target_sample_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset
    ADD CONSTRAINT callset_target_sample_id_fkey FOREIGN KEY (target_sample_id) REFERENCES sample(id);


--
-- Name: callset_to_db_array_association_callset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset_to_db_array_association
    ADD CONSTRAINT callset_to_db_array_association_callset_id_fkey FOREIGN KEY (callset_id) REFERENCES callset(id);


--
-- Name: callset_to_db_array_association_db_array_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset_to_db_array_association
    ADD CONSTRAINT callset_to_db_array_association_db_array_id_fkey FOREIGN KEY (db_array_id) REFERENCES db_array(id);


--
-- Name: callset_variant_set_callset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset_variant_set
    ADD CONSTRAINT callset_variant_set_callset_id_fkey FOREIGN KEY (callset_id) REFERENCES callset(id);


--
-- Name: callset_variant_set_variant_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY callset_variant_set
    ADD CONSTRAINT callset_variant_set_variant_set_id_fkey FOREIGN KEY (variant_set_id) REFERENCES variant_set(id);


--
-- Name: db_array_field_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY db_array
    ADD CONSTRAINT db_array_field_set_id_fkey FOREIGN KEY (field_set_id) REFERENCES field_set(id);


--
-- Name: db_array_reference_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY db_array
    ADD CONSTRAINT db_array_reference_set_id_fkey FOREIGN KEY (reference_set_id) REFERENCES reference_set(id);


--
-- Name: db_array_workspace_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY db_array
    ADD CONSTRAINT db_array_workspace_id_fkey FOREIGN KEY (workspace_id) REFERENCES workspace(id);


--
-- Name: field_field_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY field
    ADD CONSTRAINT field_field_set_id_fkey FOREIGN KEY (field_set_id) REFERENCES field_set(id);


--
-- Name: reference_reference_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference
    ADD CONSTRAINT reference_reference_set_id_fkey FOREIGN KEY (reference_set_id) REFERENCES reference_set(id);


--
-- Name: reference_set_source_accession_reference_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_set_source_accession
    ADD CONSTRAINT reference_set_source_accession_reference_set_id_fkey FOREIGN KEY (reference_set_id) REFERENCES reference_set(id);


--
-- Name: reference_set_source_accession_source_accession_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_set_source_accession
    ADD CONSTRAINT reference_set_source_accession_source_accession_id_fkey FOREIGN KEY (source_accession_id) REFERENCES source_accession(id);


--
-- Name: reference_source_accession_reference_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_source_accession
    ADD CONSTRAINT reference_source_accession_reference_id_fkey FOREIGN KEY (reference_id) REFERENCES reference(id);


--
-- Name: reference_source_accession_source_accession_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY reference_source_accession
    ADD CONSTRAINT reference_source_accession_source_accession_id_fkey FOREIGN KEY (source_accession_id) REFERENCES source_accession(id);


--
-- Name: sample_individual_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY sample
    ADD CONSTRAINT sample_individual_id_fkey FOREIGN KEY (individual_id) REFERENCES individual(id);


--
-- Name: variant_set_reference_set_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: karthikg_genomicsdb
--

ALTER TABLE ONLY variant_set
    ADD CONSTRAINT variant_set_reference_set_id_fkey FOREIGN KEY (reference_set_id) REFERENCES reference_set(id);


--
-- Name: public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- PostgreSQL database dump complete
--

