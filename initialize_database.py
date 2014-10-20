from immuno.ui import User, db, app, user_manager

db.create_all()
with app.app_context():
    admin = User(email='admin@example.com', active=True,
        password=user_manager.hash_password('1234'))
    db.session.add(admin)
    db.session.commit()
