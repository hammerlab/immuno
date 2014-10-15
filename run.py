from os.path import isdir
from immuno.ui import app

if __name__ == '__main__':
    app.debug = app.config['DEBUG']
    assert isdir(app.config['VARIANT_PATH']), \
        'Variant path %s must be a directory' % app.config['VARIANT_PATH']
    assert isdir(app.config['HLA_PATH']), \
        'HLA path %s must be a directory' % app.config['HLA_PATH']
    app.run(use_reloader=app.config['USE_RELOADER'],
            port=app.config['PORT'],
            host='0.0.0.0')
