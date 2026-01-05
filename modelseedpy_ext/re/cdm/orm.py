import json
from sqlalchemy import create_engine, Column, Integer, String, UniqueConstraint, func, insert
from sqlalchemy.dialects.sqlite import BLOB
from sqlalchemy.dialects.mysql import JSON
from sqlalchemy.orm import declarative_base, Session


Base = declarative_base()


class PanGenomeFeature(Base):

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


class GenomeFeature(Base):
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
