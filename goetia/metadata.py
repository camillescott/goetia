import datetime
import os


GOETIA_DIR = os.path.abspath(os.path.join(__file__, os.pardir))
DATA_DIR =  os.path.join(GOETIA_DIR, 'data')
TEMPLATE_DIR = os.path.join(GOETIA_DIR, 'templates')
CUR_TIME = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")

