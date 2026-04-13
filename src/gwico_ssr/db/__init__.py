"""Database engine creation and session management for GWICO-SSR."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Generator

from sqlalchemy import create_engine, event, text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session, sessionmaker

from gwico_ssr.models.schema import Base


def _set_sqlite_pragmas(dbapi_conn: object, connection_record: object) -> None:
    """Enable WAL mode and foreign keys for SQLite connections."""
    cursor = dbapi_conn.cursor()  # type: ignore[attr-defined]
    cursor.execute("PRAGMA journal_mode=WAL")
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


def get_engine(url: str, echo: bool = False) -> Engine:
    """Create a SQLAlchemy engine from a database URL.

    For SQLite URLs, automatically enables WAL mode and foreign keys.
    """
    engine = create_engine(url, echo=echo)

    if url.startswith("sqlite"):
        event.listen(engine, "connect", _set_sqlite_pragmas)

    return engine


def create_tables(engine: Engine) -> None:
    """Create all tables defined in the ORM schema."""
    Base.metadata.create_all(engine)


def drop_tables(engine: Engine) -> None:
    """Drop all tables defined in the ORM schema."""
    Base.metadata.drop_all(engine)


def get_session_factory(engine: Engine) -> sessionmaker[Session]:
    """Create a session factory bound to the given engine."""
    return sessionmaker(bind=engine, expire_on_commit=False)


@contextmanager
def session_scope(session_factory: sessionmaker[Session]) -> Generator[Session, None, None]:
    """Provide a transactional scope around a series of operations.

    Commits on clean exit, rolls back on exception.
    """
    session = session_factory()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()
