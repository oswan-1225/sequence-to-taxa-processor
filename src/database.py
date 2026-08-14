import sqlite3
import os

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