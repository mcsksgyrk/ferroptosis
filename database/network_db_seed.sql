DROP TABLE IF EXISTS `node`;
DROP TABLE IF EXISTS `edge`;
DROP TABLE IF EXISTS `aliases`;

PRAGMA foreign_keys = ON ;

CREATE TABLE `node` (
	`id`	INTEGER PRIMARY KEY,
	`name`	TEXT NOT NULL,
	`alt_accession`	TEXT,
	`tax_id`	INTEGER NOT NULL,
	`pathways`	TEXT,
	`aliases` TEXT,
	`topology` TEXT
);
CREATE TABLE `edge` (
	`id`	INTEGER PRIMARY KEY,
	`interactor_a_node_id`	INTEGER NOT NULL,
	`interactor_b_node_id`	INTEGER NOT NULL,
	`interactor_a_node_name`	TEXT NOT NULL,
	`interactor_b_node_name`	TEXT NOT NULL,
	`interaction_detection_method`	TEXT,
	`first_author`	TEXT,
	`publication_ids`	TEXT NOT NULL,
	`interaction_types`	TEXT,
	`source_db`	TEXT NOT NULL,
	`interaction_identifiers`	TEXT,
	`confidence_scores`	TEXT,
	`layer` INTEGER NOT NULL,
	FOREIGN KEY(`interactor_a_node_id`) REFERENCES node ( `id` ) ON UPDATE NO ACTION ON DELETE CASCADE,
	FOREIGN KEY(`interactor_b_node_id`) REFERENCES node ( `id` ) ON UPDATE NO ACTION ON DELETE CASCADE
);
