import duckdb
import time
import functools
import datetime
import json
from pathlib import Path


class PipelineLogger:
    """
    A class to log execution times and parameters of functions
    in a pipeline, using the @ decorator syntax.
    """

    def __init__(self, in_memory=True, db_path="pipeline_log.db"):
        """
        Initializes the DuckDB connection and the log table.

        :param in_memory: If True, uses an in-memory database (fast, volatile).
                          If False, uses a file-based database (memory-efficient, persistent).
        :param db_path: The path to the database file if in_memory is False.
        """
        db_connection = ":memory:" if in_memory else db_path
        mode = "in-memory" if in_memory else f"file-based ({db_path})"

        print(f"Initializing logger in {mode} mode...")
        self.con = duckdb.connect(database=db_connection)

        # The 'IF NOT EXISTS' clause prevents errors if the table already exists in a file-based DB
        self.con.execute("""
                         CREATE TABLE IF NOT EXISTS log_pipeline
                         (
                             nome_funzione
                             VARCHAR,
                             parametri
                             VARCHAR,
                             timestamp_inizio
                             TIMESTAMP,
                             timestamp_fine
                             TIMESTAMP,
                             durata_secondi
                             DOUBLE
                         )
                         """)

    def log(self, func):
        """
        Decorator to measure and log the execution and parameters of a function.
        """

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time_ts = time.time()

            # Capture and format parameters into a JSON string
            # The 'default=str' helps handle non-serializable data types like datetimes
            try:
                params_str = json.dumps({'args': args, 'kwargs': kwargs}, default=str)
            except TypeError:
                params_str = "Parameters not serializable to JSON"

            # Execute the original function
            result = func(*args, **kwargs)

            end_time_ts = time.time()
            duration = end_time_ts - start_time_ts

            start_time_dt = datetime.datetime.fromtimestamp(start_time_ts)
            end_time_dt = datetime.datetime.fromtimestamp(end_time_ts)

            # Insert data, including parameters, into the DuckDB table
            self.con.execute(
                "INSERT INTO log_pipeline VALUES (?, ?, ?, ?, ?)",
                (func.__name__, params_str, start_time_dt, end_time_dt, duration)
            )
            print(f"LOG: Execution of '{func.__name__}' completed in {duration:.4f} seconds.")
            return result

        return wrapper

    def save_to_csv(self, filepath="pipeline_log.csv"):
        """
        Saves the content of the log table to a CSV file.
        """
        print(f"\nSaving log to file: {filepath}...")
        try:
            # Use DuckDB's COPY command to export data to CSV
            # Use forward slashes for DuckDB compatibility on Windows
            safe_path = str(Path(filepath)).replace('\\', '/')
            self.con.execute(f"COPY log_pipeline TO '{safe_path}' (HEADER, DELIMITER ',')")
            print("Save completed successfully.")
        except Exception as e:
            print(f"Error while saving the CSV file: {e}")

    def close(self):
        """
        Closes the DuckDB connection.
        """
        if self.con:
            self.con.close()
            print("DuckDB connection closed.")