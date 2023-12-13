use std::collections::BTreeMap;
use std::collections::HashSet;
use std::env;
use std::fs::read_to_string;
use std::vec::Vec; // Instead of HashMap, so that the keys are sorted, so that testing is easier.

struct Graph {
    nodes: BTreeMap<String, HashSet<String>>,
}

fn extract_cycle(node_stack: &Vec<String>) -> String {
    let mut result = node_stack.first().unwrap().clone();
    let num_elems = node_stack.len();
    for elem in node_stack.iter().take(num_elems - 1).skip(1) {
        result.push(elem.chars().last().unwrap());
    }
    result
}

impl Graph {
    fn get_cyclic_superstrings(&mut self) -> Vec<String> {
        let mut result: Vec<String> = vec![];
        // Copying to iterate over graph node copies, otherwise graph modifies
        // itself ('pruning' operation) during iteration
        let mut all_nodes: Vec<String> = vec![];
        for focal_node in self.nodes.keys() {
            all_nodes.push(focal_node.clone())
        }
        for focal_node in all_nodes {
            let mut node_stack: Vec<String> = vec![focal_node.clone()];
            let mut cycle_stack: Vec<String> = vec![];
            let mut visited_nodes: HashSet<String> = HashSet::new();
            loop {
                let cur_node = node_stack.last().unwrap().clone();
                if visited_nodes.contains(&cur_node) {
                    if (cur_node == focal_node) {
                        // Found a 'simple cycle'
                        cycle_stack.push(cur_node.clone());
                        result.push(extract_cycle(&cycle_stack));
                        // println!("Cycle stack: {:?}", node_stack);
                        // println!{"Pre pruning:"};
                        // self.print_contents();
                        self.prune(&cycle_stack);
                        // println!{"Post pruning:"};
                        // self.print_contents();
                        break;
                    }
                    break;
                }
                let mut added_edges = 0;
                let neighbours = self.nodes.get(&cur_node);
                for neighbour in neighbours.unwrap() {
                    let num_edges_of_neighbour = self.nodes.get(neighbour).unwrap().len();
                    if num_edges_of_neighbour == 0 {
                        continue;
                    }
                    node_stack.push(neighbour.clone());
                    added_edges += 1;
                }
                if added_edges == 0 {
                    break;
                }
                visited_nodes.insert(cur_node.clone());
                cycle_stack.push(cur_node.clone());
            }
        }
        result
    }

    fn prune(&mut self, node_stack: &Vec<String>) {
        for i in 0..node_stack.len() - 1 {
            let source_node = node_stack.get(i).unwrap();
            let neighbours = self.nodes.get_mut(source_node).unwrap();
            neighbours.remove(node_stack.get(i + 1).unwrap());
        }
    }

    /// This was a previous attempt at solving the problem
    fn _get_cyclic_superstrings_2(&self) -> Vec<String> {
        let mut result: Vec<String> = vec![];
        let mut cur_nodes_to_visit: Vec<String> = vec![];
        let mut all_visited_nodes: HashSet<String> = HashSet::new();
        for focal_node in self.nodes.keys() {
            if all_visited_nodes.contains(focal_node) {
                continue;
            };
            cur_nodes_to_visit.push(focal_node.clone());
            let mut cur_visited_nodes: HashSet<String> = HashSet::new();
            let mut cyclic_superstring = String::new();
            while cur_nodes_to_visit.len() > 0 {
                let visited_node = cur_nodes_to_visit.pop().unwrap();
                all_visited_nodes.insert(visited_node.clone());
                cur_visited_nodes.insert(visited_node.clone());
                if cyclic_superstring.len() == 0 {
                    cyclic_superstring.push_str(&visited_node);
                } else {
                    cyclic_superstring.push(visited_node.chars().last().unwrap())
                };
                let neighbours = self.nodes.get(&visited_node);
                for neighbour in neighbours.unwrap() {
                    if cur_visited_nodes.contains(neighbour) {
                        result.push(cyclic_superstring.clone());
                        cur_visited_nodes.clear();
                    };
                    if all_visited_nodes.contains(neighbour) {
                        continue;
                    };
                    let num_edges_of_neighbour = self.nodes.get(neighbour).unwrap().len();
                    if num_edges_of_neighbour > 0 {
                        cur_nodes_to_visit.push(neighbour.clone());
                    }
                }
            }
        }
        result
    }

    fn shave_cyclic_superstrings(&self, cyclic_superstrings: &mut Vec<String>, k_size: usize) {
        for cyclic_superstring in cyclic_superstrings.iter_mut() {
            let new_string: String = cyclic_superstring.chars().skip(k_size - 2).collect();
            *cyclic_superstring = new_string;
        }
    }
    fn filter_cyclic_superstrings(
        &self,
        cyclic_superstrings: &Vec<String>,
        min_seq_size: usize,
    ) -> Vec<String> {
        let mut result = vec![];
        for cyclic_superstring in cyclic_superstrings {
            if cyclic_superstring.len() >= min_seq_size {
                result.push(cyclic_superstring.clone());
            }
        }
        result
    }

    /// Graph construction: add two nodes in de Bruijn graph for given k+1-mer
    fn add_nodes_from_kplus1_mer(&mut self, kplus1_mer: &str) {
        let str_size: usize = kplus1_mer.len();
        let left_kmer: String = kplus1_mer.chars().take(str_size - 1).collect();
        let right_kmer: String = kplus1_mer.chars().skip(1).collect();
        let left_kmer_match = self.nodes.get_mut(&left_kmer);
        match left_kmer_match {
            Some(value) => {
                value.insert(right_kmer.clone());
            }
            None => {
                self.nodes
                    .insert(left_kmer, HashSet::from([right_kmer.clone()]));
            }
        }
        if !self.nodes.contains_key(&right_kmer) {
            self.nodes.insert(right_kmer, HashSet::new());
        }
    }

    /// Graph construction: builds the de Bruijn graph for input set of `sequences`
    fn add_kplus1_mers_from_sequences(
        &mut self,
        sequences: &Vec<String>,
        k_size: usize,
        seq_size: usize,
    ) {
        for sequence in sequences {
            let revcomp = revcomp(sequence);
            for i in 0..seq_size - k_size + 1 {
                self.add_nodes_from_kplus1_mer(&sequence[i..i + k_size]);
                self.add_nodes_from_kplus1_mer(&revcomp[i..i + k_size]);
            }
        }
    }
    /// Graph introspection
    fn print_contents(&self) {
        for (node, neighbours) in &self.nodes {
            println!("Node: {node}, neighbours: {:?}", neighbours);
        }
    }
}

fn load_sequences(filename: &str) -> Vec<String> {
    read_to_string(filename)
        .unwrap()
        .lines()
        .map(String::from)
        .collect()
}

fn revcomp(input_string: &str) -> String {
    let mut result = String::new();
    for s in input_string.to_uppercase().chars().rev() {
        match s {
            'A' => result.push('T'),
            'C' => result.push('G'),
            'G' => result.push('C'),
            'T' => result.push('A'),
            _ => panic!("Unrecognised DNA character: {s}"),
        };
    }
    result
}

fn is_cyclic_permutation(query_str: &str, target_str: &str) -> bool {
    let mut doubled_string: String = target_str.chars().collect();
    doubled_string.push_str(target_str.clone());
    return doubled_string.contains(query_str);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let sequences = load_sequences(&args[1]);
    let seq_size: usize = sequences[0].len();
    let mut template_superstring = String::new();
    // println!("Input sequences: {:?}", sequences);
    for kplus1_mer_size in (3..seq_size + 1).rev() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([]),
        };
        graph.add_kplus1_mers_from_sequences(&sequences, kplus1_mer_size, seq_size);
        // graph.print_contents();
        let mut cyclic_superstrings = graph.get_cyclic_superstrings();
        graph.shave_cyclic_superstrings(&mut cyclic_superstrings, kplus1_mer_size);
        // Rationale for filtering: a superstring must be able to contain all the reads
        cyclic_superstrings = graph.filter_cyclic_superstrings(&cyclic_superstrings, seq_size);
        // println!("k+1 size: {kplus1_mer_size}; Num cyclic superstrings {:?}", cyclic_superstrings.len());
        if cyclic_superstrings.len() == 2 {
            println!("Two cyclic superstrings found, at k+1-mer size {kplus1_mer_size}!");
            println!("These are: {:?}", cyclic_superstrings);
            if template_superstring.len() == 0 {
                template_superstring = cyclic_superstrings[0].clone();
            } else {
                let identical_superstrings_found: Vec<bool> = cyclic_superstrings
                    .iter()
                    .map(|q| is_cyclic_permutation(q, &template_superstring))
                    .collect();
                if !identical_superstrings_found.contains(&true) {
                    println!("Warning: cyclic superstrings for k+1-mer size of {kplus1_mer_size} found not identical to template");
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_node_to_empty_graph() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([]),
        };
        graph.add_nodes_from_kplus1_mer(&String::from("ATGG"));
        assert_eq!(graph.nodes.len(), 2);
        let result = graph.nodes.get("ATG");
        assert_eq!(result.unwrap().len(), 1);
        assert!(result.unwrap().contains(&String::from("TGG")));
    }

    #[test]
    fn test_add_existing_kplus1_mer_to_graph_does_nothing() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (String::from("ATG"), HashSet::from([String::from("TGC")])),
                (String::from("TGC"), HashSet::new()),
            ]),
        };
        graph.add_nodes_from_kplus1_mer(&String::from("ATGC"));
        assert_eq!(graph.nodes.len(), 2);
    }

    #[test]
    fn test_add_kplus1_mer_to_graph_with_existing_node() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (String::from("ATG"), HashSet::from([String::from("TGC")])),
                (String::from("TGC"), HashSet::new()),
            ]),
        };
        graph.add_nodes_from_kplus1_mer(&String::from("ATGA"));
        assert_eq!(graph.nodes.len(), 3);
        let result = graph.nodes.get("ATG");
        assert_eq!(result.unwrap().len(), 2);
        assert!(result.unwrap().contains(&String::from("TGC")));
        assert!(result.unwrap().contains(&String::from("TGA")));
    }

    #[test]
    fn test_no_directed_cycle() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (String::from("ATG"), HashSet::from([String::from("TGC")])),
                (String::from("TGC"), HashSet::new()),
            ]),
        };
        assert_eq!(0, graph.get_cyclic_superstrings().len());
    }

    #[test]
    fn test_single_directed_cycle() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (String::from("ATG"), HashSet::from([String::from("TGA")])),
                (String::from("TGA"), HashSet::from([String::from("GAT")])),
                (String::from("GAT"), HashSet::from([String::from("ATG")])),
            ]),
        };
        assert_eq!(vec![String::from("ATGAT")], graph.get_cyclic_superstrings());
    }

    #[test]
    fn test_single_directed_cycle_with_empty_branchout() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (String::from("ATG"), HashSet::from([String::from("TGA")])),
                (
                    String::from("TGA"),
                    HashSet::from([String::from("GAT"), String::from("GAC")]),
                ),
                (String::from("GAC"), HashSet::new()),
                (String::from("GAT"), HashSet::from([String::from("ATG")])),
            ]),
        };
        assert_eq!(vec![String::from("ATGAT")], graph.get_cyclic_superstrings());
    }

    #[test]
    fn test_two_directed_cycles() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (String::from("ATG"), HashSet::from([String::from("TGA")])),
                (String::from("TGA"), HashSet::from([String::from("GAT")])),
                (String::from("GAT"), HashSet::from([String::from("ATG")])),
                (String::from("TAA"), HashSet::from([String::from("AAA")])),
                (String::from("AAA"), HashSet::from([String::from("AAT")])),
                (String::from("AAT"), HashSet::from([String::from("ATA")])),
                (String::from("ATA"), HashSet::from([String::from("TAA")])),
            ]),
        };
        assert_eq!(
            vec![String::from("AAATAA"), String::from("ATGAT")],
            graph.get_cyclic_superstrings()
        );
    }

    #[test]
    fn test_two_directed_cycles_with_one_shared_node() {
        let mut graph: Graph = Graph {
            nodes: BTreeMap::from([
                (
                    String::from("ATG"),
                    HashSet::from([String::from("TGA"), String::from("TGC")]),
                ),
                (String::from("TGA"), HashSet::from([String::from("GAT")])),
                (String::from("GAT"), HashSet::from([String::from("ATG")])),
                (String::from("TGC"), HashSet::from([String::from("GCA")])),
                (String::from("GCA"), HashSet::from([String::from("CAT")])),
                (String::from("CAT"), HashSet::from([String::from("ATG")])),
            ]),
        };
        // graph.print_contents();
        assert_eq!(
            vec![String::from("ATGCAT"), String::from("GATGA")],
            graph.get_cyclic_superstrings()
        );
    }

    #[test]
    fn test_revcomp() {
        let input_str = String::from("TGTAA");
        assert_eq!(revcomp(&input_str), String::from("TTACA"));
    }

    #[test]
    fn test_cyclic_permutations() {
        assert!(is_cyclic_permutation(
            &String::from("TTACA"),
            &String::from("ACATT")
        ));
        assert!(!is_cyclic_permutation(
            &String::from("TTAA"),
            &String::from("TTTT")
        ));
    }
}
