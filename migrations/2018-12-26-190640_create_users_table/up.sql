-- This creates the users table with all the required information for session management

CREATE EXTENSION "uuid-ossp";

CREATE TABLE users (
  id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
  username VARCHAR(35), -- Username of the user
  password CHAR(70), -- 64 characters for SHA-256, 6 for the salt
  email VARCHAR(65), -- Email of the user
  reg_time TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP, -- Time the account was created
  last_online TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP -- Time of the last connection of the account
);
