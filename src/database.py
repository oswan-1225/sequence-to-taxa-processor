import sqlite3
import os
import pandas as pd

def create_database(db_path: str) -> None:
    """
    Create an SQLite database and its tables, if they do not already exist.
    
    Parameters:
        db_path (str): The path to the SQLite database file.
    """
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS samples (
            sample_id TEXT PRIMARY KEY,
            source TEXT,
            notes TEXT      
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS classifications (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            read_id TEXT,
            best_match TEXT,
            confidence REAL,
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
        )
        """)

    conn.commit()
    conn.close()

def insert_sample_results(db_path: str, sample_id: str, source: str, results: list[dict]) -> None:
    """
    Insert a sample's classifciation results into the database.
    Creates the sample rwo if it doesn't already exist.
    
    Parameters:
        db_path (str): Path to the SQLite database file.
        sample_id (str): Unique identifier for the sample.
        source (str): Source of the sample (e.g., "SRA", "local").
        results (list[dict]): List of {'read_id', "best_match", "confidence"} dicts, matching the shape classify_reads.py already builds.
        """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("PRAGMA foreign_keys = ON")

    cursor.execute("""INSERT OR IGNORE INTO samples (sample_id, source) VALUES (?, ?)""", (sample_id, source))

    classification_rows = [(sample_id, result['read_id'], result['best_match'], result['confidence']) for result in results]
    cursor.executemany("""INSERT INTO classifications (sample_id, read_id, best_match, confidence) VALUES (?, ?, ?, ?)""", classification_rows)
    conn.commit()
    conn.close()

def get_abundance(db_path: str) -> pd.DataFrame:
    """
    Retrieve the abundance of each taxon across all samples in the database.
    
    Parameters:
        db_path (str): Path to the SQLite database file.
    
    Returns:
        pd.DataFrame: columns ['sample_id', 'best_match', 'count', 'total', 'percent'] - one row
                      per (sample, species), with 'total' the sample's classified read count and
                      'percent' the species' share of that total.
    """
    conn = sqlite3.connect(db_path)

    query = """
        SELECT sample_id, best_match, COUNT(*) AS count
        FROM classifications
        WHERE best_match IS NOT NULL
        GROUP BY sample_id, best_match
    """
    df = pd.read_sql_query(query, conn)

    df['total'] = df.groupby('sample_id')['count'].transform('sum')
    df['percent'] = df['count'] / df['total'] * 100

    conn.close()
    return df

def get_classification_totals(db_path: str) -> pd.DataFrame:
    """
    Per-sample read counts, including reads with no k-mer match to any
    reference species (best_match IS NULL) - the rows get_abundance()
    excludes.

    Parameters:
        db_path (str): Path to the SQLite database file.

    Returns:
        pd.DataFrame: columns ['sample_id', 'total_reads', 'classified_reads',
                      'unclassified_reads', 'unclassified_percent'].
    """
    conn = sqlite3.connect(db_path)

    query = """
        SELECT sample_id,
               COUNT(*) AS total_reads,
               SUM(CASE WHEN best_match IS NOT NULL THEN 1 ELSE 0 END) AS classified_reads
        FROM classifications
        GROUP BY sample_id
    """
    df = pd.read_sql_query(query, conn)
    conn.close()

    df['unclassified_reads'] = df['total_reads'] - df['classified_reads']
    df['unclassified_percent'] = df['unclassified_reads'] / df['total_reads'] * 100
    return df