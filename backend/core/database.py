from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from backend.core.config import get_settings

settings = get_settings()

# Create the SQLAlchemy engine
# connect_args={"check_same_thread": False} is needed for SQLite
# to allow multiple threads to interact with the same connection.
engine = create_engine(
    settings.DATABASE_URL, connect_args={"check_same_thread": False}
)

# Create a SessionLocal class
# Each instance of SessionLocal will be a database session.
# The expire_on_commit=False will prevent objects from being expired
# after commit, so you can access their attributes after committing.
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Base class for our models to inherit from
Base = declarative_base()

# Dependency to get a database session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
