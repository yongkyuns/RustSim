//! Subsystem navigation for hierarchical graph traversal.

use rustsim_types::SimulationGraph;

/// Navigates through nested subsystems in a hierarchical graph
#[derive(Debug, Clone, Default)]
pub struct SubsystemNavigator {
    /// Stack of subsystem node IDs (empty = root level)
    path: Vec<String>,
}

impl SubsystemNavigator {
    /// Create a new navigator at the root level
    pub fn new() -> Self {
        Self::default()
    }

    /// Get the current graph at the current level of navigation
    pub fn current_graph<'a>(&self, root: &'a SimulationGraph) -> Result<&'a SimulationGraph, String> {
        self.traverse_to_current(root)
    }

    /// Get the current graph mutably at the current level of navigation
    pub fn current_graph_mut<'a>(
        &self,
        root: &'a mut SimulationGraph,
    ) -> Result<&'a mut SimulationGraph, String> {
        self.traverse_to_current_mut(root)
    }

    /// Enter a subsystem by its node ID
    pub fn enter_subsystem(&mut self, subsystem_id: &str) -> Result<(), String> {
        // Validation would happen in the caller (checking if the node exists and is a subsystem)
        self.path.push(subsystem_id.to_string());
        Ok(())
    }

    /// Exit the current subsystem and go up one level
    pub fn exit_subsystem(&mut self) -> Result<(), String> {
        if self.path.is_empty() {
            return Err("Already at root level".to_string());
        }
        self.path.pop();
        Ok(())
    }

    /// Go back to the root level
    pub fn go_to_root(&mut self) {
        self.path.clear();
    }

    /// Check if currently at root level
    pub fn is_at_root(&self) -> bool {
        self.path.is_empty()
    }

    /// Get the current depth (0 = root)
    pub fn current_depth(&self) -> usize {
        self.path.len()
    }

    /// Get breadcrumb path as (id, name) pairs for UI display
    pub fn breadcrumb_path(&self, root: &SimulationGraph) -> Vec<(String, String)> {
        let mut breadcrumbs = Vec::new();
        let current_graph = root;

        for subsystem_id in &self.path {
            // Get the node to extract its name
            if let Some(node) = current_graph.get_node(subsystem_id) {
                breadcrumbs.push((subsystem_id.clone(), node.name.clone()));

                // Move to the nested graph for the next iteration
                if let Some(_subsystem_graph) = &node.graph {
                    // Convert SubsystemGraph to SimulationGraph for traversal
                    // This is a workaround since SubsystemGraph stores nodes as Vec instead of HashMap
                    // In practice, we'd need to convert or have a better structure
                    // For now, we'll stop traversing deeper for breadcrumb name resolution
                    break;
                }
            } else {
                // If we can't find the node, use the ID as the name
                breadcrumbs.push((subsystem_id.clone(), subsystem_id.clone()));
                break;
            }
        }

        breadcrumbs
    }

    /// Get the full path as a slice
    pub fn path(&self) -> &[String] {
        &self.path
    }

    /// Helper to traverse to the current level (immutable)
    fn traverse_to_current<'a>(&self, root: &'a SimulationGraph) -> Result<&'a SimulationGraph, String> {
        if self.path.is_empty() {
            return Ok(root);
        }

        let current_graph = root;

        for subsystem_id in &self.path {
            let node = current_graph
                .get_node(subsystem_id)
                .ok_or_else(|| format!("Subsystem node '{}' not found", subsystem_id))?;

            if node.block_type != "Subsystem" {
                return Err(format!("Node '{}' is not a Subsystem", subsystem_id));
            }

            // SubsystemGraph needs to be converted to SimulationGraph
            // Since NodeInstance has Option<Box<SubsystemGraph>>, we need to handle the conversion
            // For now, this will be a limitation - we'll return an error
            return Err("Nested subsystem navigation not yet fully implemented".to_string());
        }

        Ok(current_graph)
    }

    /// Helper to traverse to the current level (mutable)
    fn traverse_to_current_mut<'a>(
        &self,
        root: &'a mut SimulationGraph,
    ) -> Result<&'a mut SimulationGraph, String> {
        if self.path.is_empty() {
            return Ok(root);
        }

        let current_graph = root;

        for subsystem_id in &self.path {
            let node = current_graph
                .get_node_mut(subsystem_id)
                .ok_or_else(|| format!("Subsystem node '{}' not found", subsystem_id))?;

            if node.block_type != "Subsystem" {
                return Err(format!("Node '{}' is not a Subsystem", subsystem_id));
            }

            // SubsystemGraph needs to be converted to SimulationGraph
            // This is a limitation of the current design
            return Err("Nested subsystem navigation not yet fully implemented".to_string());
        }

        Ok(current_graph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsim_types::{NodeInstance, Position};

    #[test]
    fn test_navigator_basic() {
        let nav = SubsystemNavigator::new();
        assert!(nav.is_at_root());
        assert_eq!(nav.current_depth(), 0);
        assert_eq!(nav.path().len(), 0);
    }

    #[test]
    fn test_enter_exit_subsystem() {
        let mut nav = SubsystemNavigator::new();

        assert!(nav.enter_subsystem("sub1").is_ok());
        assert!(!nav.is_at_root());
        assert_eq!(nav.current_depth(), 1);

        assert!(nav.enter_subsystem("sub2").is_ok());
        assert_eq!(nav.current_depth(), 2);

        assert!(nav.exit_subsystem().is_ok());
        assert_eq!(nav.current_depth(), 1);

        assert!(nav.exit_subsystem().is_ok());
        assert!(nav.is_at_root());

        // Try to exit when already at root
        assert!(nav.exit_subsystem().is_err());
    }

    #[test]
    fn test_go_to_root() {
        let mut nav = SubsystemNavigator::new();

        nav.enter_subsystem("sub1").unwrap();
        nav.enter_subsystem("sub2").unwrap();
        nav.enter_subsystem("sub3").unwrap();
        assert_eq!(nav.current_depth(), 3);

        nav.go_to_root();
        assert!(nav.is_at_root());
        assert_eq!(nav.current_depth(), 0);
    }

    #[test]
    fn test_breadcrumb_path_at_root() {
        let nav = SubsystemNavigator::new();
        let graph = SimulationGraph::new();

        let breadcrumbs = nav.breadcrumb_path(&graph);
        assert_eq!(breadcrumbs.len(), 0);
    }

    #[test]
    fn test_breadcrumb_path_with_subsystems() {
        let mut nav = SubsystemNavigator::new();
        let mut graph = SimulationGraph::new();

        // Add a subsystem node
        let mut node = NodeInstance::new(
            "sub1".to_string(),
            "Subsystem".to_string(),
            Position::zero(),
        );
        node.name = "My Subsystem".to_string();
        graph.add_node(node);

        nav.enter_subsystem("sub1").unwrap();

        let breadcrumbs = nav.breadcrumb_path(&graph);
        assert_eq!(breadcrumbs.len(), 1);
        assert_eq!(breadcrumbs[0].0, "sub1");
        assert_eq!(breadcrumbs[0].1, "My Subsystem");
    }
}
