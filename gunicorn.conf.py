"""Production server settings. gunicorn reads this file from the working directory."""
import os

bind = f"0.0.0.0:{os.environ.get('PORT', '8000')}"
workers = 1               # one copy of the app fits the free tier's 512 MB
threads = 4               # pages and quick requests still answer while a search runs
timeout = 120             # an alignment search takes seconds locally, longer on a 0.1 CPU share
limit_request_line = 0    # a link carrying a 720-point skyline is about 10 kB; the default cap is 4 kB
accesslog = '-'           # one line per request in the host's log
