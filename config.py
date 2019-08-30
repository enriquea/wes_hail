import os

# The absolute path of this project
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# The logs files path
LOGS_PATH = os.path.join(ROOT_DIR, 'logs')

# Hail log file
HAIL_LOG_PATH = os.path.join(ROOT_DIR, 'logs/hail.log')

# Hail temporal directory path (network visible)
HAIL_TMP_DIR = '/mnt/nfs/mdatanode/hail-temp'
