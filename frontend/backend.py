from flask import Flask, jsonify, send_from_directory
import sqlite3

app = Flask(__name__)


def get_data():
    conn = sqlite3.connect('./omnipath2.db')
    cursor = conn.cursor()

    # Fetch nodes
    cursor.execute("SELECT id, name, pathways, function FROM node")
    nodes = [{"id": row[0], "name": row[1], "pathways": row[2], "function": row[3]} for row in cursor.fetchall()]

    # Fetch edges
    cursor.execute("SELECT interactor_a_node_id, interactor_b_node_id, layer FROM edge")
    edges = [{"source": row[0], "target": row[1], "layer": row[2]} for row in cursor.fetchall()]

    conn.close()
    return {"nodes": nodes, "edges": edges}


@app.route('/')
def serve_index():
    return send_from_directory('templates', 'index.html')


@app.route('/static/<path:path>')
def serve_static(path):
    return send_from_directory('static', path)


@app.route('/graph-data')
def graph_data():
    return jsonify(get_data())


if __name__ == '__main__':
    app.run(debug=True)
