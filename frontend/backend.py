from flask import Flask, jsonify, send_from_directory
import sqlite3

app = Flask(__name__)


def get_data():
    conn = sqlite3.connect('graph.db')
    cursor = conn.cursor()

    # Fetch nodes
    cursor.execute("SELECT id, label, 'group' FROM nodes")
    nodes = [{"id": row[0], "label": row[1], "group": row[2]} for row in cursor.fetchall()]

    # Fetch edges
    cursor.execute("SELECT source, target, weight FROM edges")
    edges = [{"source": row[0], "target": row[1], "weight": row[2]} for row in cursor.fetchall()]

    conn.close()
    return {"nodes": nodes, "edges": edges}


@app.route('/graph-data')
def graph_data():
    return jsonify(get_data())


@app.route('/')
def serve_frontend():
    return send_from_directory('.', 'index.html')


if __name__ == '__main__':
    app.run(debug=True)
