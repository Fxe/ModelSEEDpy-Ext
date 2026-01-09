import json
from sqlalchemy import create_engine, Column, Integer, Text, String, Float, Boolean, UniqueConstraint, func, insert
from sqlalchemy.dialects.sqlite import BLOB
from sqlalchemy.dialects.mysql import JSON
from sqlalchemy.orm import declarative_base, Session


Base = declarative_base()


class PhenotypeKeggModule(Base):
    __tablename__ = "phenotype_module"

    id         = Column(String(6), primary_key=True)
    label  = Column(String(500), nullable=False)
    definition  = Column(String(1000), nullable=False)
    kos = Column(String(1000), nullable=False)
    n_steps = Column(Integer, nullable=False)
    steps_in_genome = Column(Integer, nullable=False)
    steps_in_pangenome = Column(Integer, nullable=False)
    steps_in_pangenome_hybrid = Column(Integer, nullable=False)
    frac_in_genome = Column(Float, nullable=False)
    frac_in_pangenome = Column(Float, nullable=False)
    frac_in_pangenome_hybrid = Column(Float, nullable=False)
    frac_any = Column(Float, nullable=False)


class GenomeMetadata(Base):
    __tablename__ = "genome"

    id = Column(String(255), primary_key=True)
    gtdb_taxonomy = Column(String(1000))
    ncbi_taxonomy = Column(String(1000))
    n_contigs = Column(Integer)
    n_features = Column(Integer)


class GenomeANI(Base):
    __tablename__ = "genome_ani"

    genome1 = Column(String(255), primary_key=True)
    genome2 = Column(String(255), primary_key=True)

    ani = Column(Float, nullable=False)
    af1 = Column(Float, nullable=False)
    af2 = Column(Float, nullable=False)

    kind = Column(String(255), nullable=False)


class GenomeFeature(Base):
    __tablename__ = "genome_features"

    id = Column(Integer, primary_key=True)
    genome_id = Column(String(255), nullable=False)
    contig_id = Column(String(255), nullable=False)
    feature_id = Column(String(255), nullable=False)
    length = Column(Integer, nullable=False)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String(1))
    sequence = Column(Text)
    sequence_hash = Column(String(255))

    bakta_function = Column(String(255))
    rast_function = Column(String(255))
    cog = Column(String(255))
    ec = Column(String(255))
    gene_names = Column(String(255))

    go = Column(String(1000))
    ko = Column(String(1000))
    pfam = Column(String(1000))
    go = Column(String(1000))
    so = Column(String(1000))

    uniref_100 = Column(String(255))
    uniref_90 = Column(String(255))
    uniref_50 = Column(String(255))

    pangenome_cluster_id = Column(String(255))
    pangenome_is_core = Column(Boolean)

    psortb = Column(String(255))

    __table_args__ = (
        UniqueConstraint("genome_id", "contig_id", "feature_id", name="uq_feature"),
    )


class PanGenomeFeature(Base):

    __tablename__ = "pan_genome_features"

    id         = Column(Integer, primary_key=True)
    genome_id  = Column(String(255), nullable=False)
    contig_id  = Column(String(255), nullable=False)
    feature_id = Column(String(255), nullable=False)
    length     = Column(Integer, nullable=False)
    start      = Column(Integer)
    end        = Column(Integer)
    strand     = Column(String(1))
    sequence   = Column(Text)
    sequence_hash   = Column(String(255))

    cluster_id = Column(String(255))
    is_core    = Column(Boolean)

    bakta_function  = Column(String(255))
    rast_function   = Column(String(255))

    gene_names      = Column(String(255))

    cog             = Column(String(1000))
    ec              = Column(String(1000))
    ko     = Column(String(1000))
    pfam   = Column(String(1000))
    go     = Column(String(1000))
    so     = Column(String(1000))

    uniref_100    = Column(String(255))
    uniref_90     = Column(String(255))
    uniref_50     = Column(String(255))
    #annotation = Column(BLOB, nullable=False)  # JSONB stored as BLOB
    #annotation = Column(JSON, nullable=False)  # MySQL native JSON

    __table_args__ = (
        UniqueConstraint("genome_id", "contig_id", "feature_id", name="uq_feature"),
    )


BaseL = declarative_base()


class PanGenomeFeatureL(BaseL):

    __tablename__ = "pan_genome_features"

    id         = Column(Integer, primary_key=True)
    genome_id  = Column(String, nullable=False)
    contig_id  = Column(String, nullable=False)
    feature_id = Column(String, nullable=False)
    length     = Column(Integer, nullable=False)
    start      = Column(Integer)
    end        = Column(Integer)
    strand     = Column(String)
    annotation = Column(BLOB, nullable=False)  # JSONB stored as BLOB
    #annotation = Column(JSON, nullable=False)  # MySQL native JSON

    __table_args__ = (
        UniqueConstraint("genome_id", "contig_id", "feature_id", name="uq_feature"),
    )


class GenomeFeatureL(BaseL):
    __tablename__ = "genome_features"

    id         = Column(Integer, primary_key=True)
    genome_id  = Column(String, nullable=False)
    contig_id  = Column(String, nullable=False)
    feature_id = Column(String, nullable=False)
    length     = Column(Integer, nullable=False)
    start      = Column(Integer)
    end        = Column(Integer)
    strand     = Column(String)
    annotation = Column(BLOB, nullable=False)  # JSONB stored as BLOB

    __table_args__ = (
        UniqueConstraint("genome_id", "contig_id", "feature_id", name="uq_feature"),
    )


def load_tables(table_pan_genome_features):
    engine = create_engine("sqlite:///berdl_tables.db", future=True)

    Base.metadata.create_all(engine)
    rows = [
        dict(
            genome_id=r.genome_id,
            contig_id=r.contig_id,
            feature_id=r.feature_id,
            length=int(r.length),
            start=r.start,
            end=r.end,
            strand=r.strand,
            annotation=json.dumps(r.annotation, ensure_ascii=False),
        )
        for r in table_pan_genome_features
    ]

    # Use Core insert, so we can apply jsonb() function
    stmt = insert(PanGenomeFeature).values([
        {**row, "annotation": func.jsonb(row["annotation"])} for row in rows
    ])

    with engine.begin() as conn:
        conn.execute(stmt)


def load_tables_mysql(table_pan_genome_features):
    engine = create_engine(
        "mysql+pymysql://USER:PASSWORD@HOST:3306/DBNAME?charset=utf8mb4",
        future=True
    )
    Base.metadata.create_all(engine)

    with Session(engine) as ses:
        ses.add_all([
            PanGenomeFeature(
                genome_id=r.genome_id,
                contig_id=r.contig_id,
                feature_id=r.feature_id,
                length=int(r.length),
                start=r.start,
                end=r.end,
                strand=r.strand,
                annotation=r.annotation,  # pass dict directly
            )
            for r in table_pan_genome_features
        ])
        ses.commit()
