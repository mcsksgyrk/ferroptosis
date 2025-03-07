DROP TABLE IF EXISTS `node`;
DROP TABLE IF EXISTS `edge`;
DROP TABLE IF EXISTS `node_identifier`;
PRAGMA foreign_keys = ON;

CREATE TABLE `node` (
    `id` INTEGER PRIMARY KEY,
    `name` TEXT NOT NULL,         -- Primary identifier (UniProt/PubChem/KEGG ID)
    `display_name` TEXT,          -- Human-readable name
    `tax_id` INTEGER,             -- Only for biological entities, now nullable
    `type` TEXT NOT NULL DEFAULT 'protein', -- 'protein', 'small_molecule', 'rna', etc.
    `pathways` TEXT,
    `source` TEXT,
    `function` TEXT
);

CREATE TABLE `node_identifier` (
    `node_id` INTEGER NOT NULL,
    `id_type` TEXT NOT NULL,      -- 'kegg_id', 'uniprot_id', 'pubchem_id', etc.
    `id_value` TEXT NOT NULL,
    PRIMARY KEY (`node_id`, `id_type`),
    FOREIGN KEY (`node_id`) REFERENCES `node`(`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

CREATE TABLE `edge` (
    `id` INTEGER PRIMARY KEY,
    `interactor_a_node_id` INTEGER NOT NULL,
    `interactor_b_node_id` INTEGER NOT NULL,
    `interactor_a_node_name` TEXT NOT NULL,
    `interactor_b_node_name` TEXT NOT NULL,
    `layer` TEXT NOT NULL,
    `interaction_types` TEXT,
    `source_db` TEXT NOT NULL,
    FOREIGN KEY(`interactor_a_node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE,
    FOREIGN KEY(`interactor_b_node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

-- Add indexes to improve query performance
CREATE INDEX `idx_node_name` ON `node`(`name`);
CREATE INDEX `idx_node_name_tax` ON `node`(`name`, `tax_id`);
CREATE INDEX `idx_node_type` ON `node`(`type`);
CREATE INDEX `idx_node_identifier_type` ON `node_identifier`(`id_type`);
CREATE INDEX `idx_node_identifier_value` ON `node_identifier`(`id_value`);
CREATE INDEX `idx_node_identifier_node` ON `node_identifier`(`node_id`);
CREATE INDEX `idx_edge_interactors` ON `edge`(`interactor_a_node_id`, `interactor_b_node_id`);
CREATE INDEX `idx_edge_layer` ON `edge`(`layer`);
