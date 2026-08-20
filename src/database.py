import sqlite3
import os
import pandas as pd

# Every column create_database() expects on the classifications table. Used to
# detect a database written before a schema change - see create_database().
EXPECTED_CLASSIFICATION_COLUMNS = {
    "id", "sample_id", "read_id", "best_match", "confidence", "n_species_hit"
}


def create_database(db_path: str) -> None:
    """
    Create an SQLite database and its tables, if they do not already exist.

    Raises ValueError if db_path is an existing database whose classifications
    table predates a schema change. CREATE TABLE IF NOT EXISTS is a no-op on an
    existing table, so without this check an outdated database silently keeps
    its old columns and the failure surfaces later as an opaque sqlite error
    inside executemany. Results databases are regenerable, so the fix is to
    delete the file and rerun.

    Parameters:
        db_path (str): The path to the SQLite database file.
    """
    # abspath first: for a bare filename like "results.db", os.path.dirname
    # returns "" and os.makedirs("") raises FileNotFoundError.
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
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
            n_species_hit INTEGER,
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
        )
        """)

    # Runs after the CREATE, which is a no-op when the table already exists.
    # On a fresh database the statement above just created every expected
    # column, so this passes; on an outdated one it fires.
    cursor.execute("PRAGMA table_info(classifications)")
    existing_columns = {row[1] for row in cursor.fetchall()}
    missing = EXPECTED_CLASSIFICATION_COLUMNS - existing_columns
    if missing:
        conn.close()
        raise ValueError(
            f"{db_path} has an outdated classifications table, missing "
            f"{sorted(missing)}. Classification databases are regenerable "
            f"outputs, so delete it and rerun rather than migrating it."
        )

    conn.commit()
    conn.close()

def insert_sample_results(db_path: str, sample_id: str, source: str, results: list[dict]) -> None:
    """
    Insert a sample's classification results into the database.
    Creates the sample row if it doesn't already exist.

    Re-inserting an existing sample_id DELETES that sample's previous
    classification rows before writing the new ones, so running the same
    sample twice gives the same result as running it once rather than
    doubling its counts. Other samples in the database are untouched.

    Parameters:
        db_path (str): Path to the SQLite database file.
        sample_id (str): Unique identifier for the sample.
        source (str): Source of the sample (e.g., "SRA", "local").
        results (list[dict]): List of {'read_id', "best_match", "confidence",
            "n_species_hit"} dicts, matching the shape classify_reads.py
            already builds. n_species_hit is read with .get() rather than [],
            so callers that predate the column (and tests that hand-build
            result dicts) still work, storing NULL for it.
        """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("PRAGMA foreign_keys = ON")

    cursor.execute("""INSERT OR IGNORE INTO samples (sample_id, source) VALUES (?, ?)""", (sample_id, source))

    cursor.execute("DELETE FROM classifications WHERE sample_id = ?", (sample_id,))

    classification_rows = [(sample_id, result['read_id'], result['best_match'],
                            result['confidence'], result.get('n_species_hit'))
                           for result in results]
    cursor.executemany(
        """INSERT INTO classifications (sample_id, read_id, best_match, confidence, n_species_hit)
           VALUES (?, ?, ?, ?, ?)""",
        classification_rows)
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