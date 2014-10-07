from flask import Flask
app = Flask(__name__)

@app.route('/')
def index():
    return 'Immuno UI'
	
@app.route('/patients')
def patients():
    return 'Patients'

if __name__ == '__main__':
    app.run()
