import datetime
import os


BOINK_DIR = os.path.abspath(os.path.join(__file__, os.pardir))
DATA_DIR =  os.path.join(BOINK_DIR, 'data')
TEMPLATE_DIR = os.path.join(BOINK_DIR, 'templates')
CUR_TIME = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")

