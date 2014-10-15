from os.path import isdir
from immuno.ui import app

if __name__ == '__main__':
    app.debug = app.config['DEBUG']
    app.run(use_reloader=app.config['USE_RELOADER'],
            port=app.config['PORT'],
            host='0.0.0.0')
