use std::env;
use clam::prelude::*;
use random_string::generate;

fn main() {
    let args: Vec<String> = env::args().collect();
    let algorithm = &args[1];
    let seq1 = &args[2];
    let seq2= &args[3];
    if algorithm == "Needleman-Wunsch" {
        algorithms_lib::Needleman_Wunsch::align(seq1.to_string(), seq2.to_string());
    }
    else if algorithm == "Smith-Waterman" {
        algorithms_lib::Smith_Waterman::align(seq1.to_string(), seq2.to_string());
    }
    else {
        println!("Algorithm is unknown. Please run program again!");
    }
}

#[cfg(test)]
mod test {
    mod test_nw {
        mod test_nw_5 {
            use clam::Metric;
            use algorithms_lib::Needleman_Wunsch::NeedlemanWunsch;
            
            #[test]
            fn metric_nw_5to5_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test] 
            fn metric_nw_5to10_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to20_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to30_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to40_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_5to50_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_10 {
            use clam::Metric;
            use algorithms_lib::Needleman_Wunsch::NeedlemanWunsch;            
            #[test]
            fn metric_nw_10to5_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test] 
            fn metric_nw_10to10_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to20_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to30_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to40_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_10to50_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_20 {
            use clam::Metric;
            use algorithms_lib::Needleman_Wunsch::NeedlemanWunsch;
            #[test]
            fn metric_nw_20to5_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test] 
            fn metric_nw_20to10_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to20_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to30_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to40_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_20to50_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_30 {
            use clam::Metric;
            use algorithms_lib::Needleman_Wunsch::NeedlemanWunsch;
            #[test]
            fn metric_nw_30to5_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test] 
            fn metric_nw_30to10_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to20_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to30_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to40_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_30to50_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }

        mod test_nw_40 {
            use clam::Metric;
            use algorithms_lib::Needleman_Wunsch::NeedlemanWunsch;
            #[test]
            fn metric_nw_40to5_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test] 
            fn metric_nw_40to10_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to20_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to30_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to40_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }

            #[test]
            fn metric_nw_40to50_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    }

    mod test_nw_50 {
        use clam::Metric;
        use algorithms_lib::Needleman_Wunsch::NeedlemanWunsch;
        #[test]
        fn metric_nw_50to5_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "CAACC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test] 
        fn metric_nw_50to10_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TTGCTTTGAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to20_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to30_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
            let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to40_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }

        #[test]
        fn metric_nw_50to50_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: NeedlemanWunsch = NeedlemanWunsch{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    }


    mod test_sw {
        mod test_sw_5 {
            use clam::Metric;
            use algorithms_lib::Smith_Waterman::SmithWaterman;
            
            #[test]
            fn metric_sw_5to5_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test] 
            fn metric_sw_5to10_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_5to20_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_5to30_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_5to40_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_5to50_test(){
                let seq1_c: String = "AAACT".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    
        mod test_sw_10 {
            use clam::Metric;
            use algorithms_lib::Smith_Waterman::SmithWaterman;            
            #[test]
            fn metric_sw_10to5_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test] 
            fn metric_sw_10to10_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_10to20_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_10to30_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_10to40_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_10to50_test(){
                let seq1_c: String = "CAGAATATTA".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    
        mod test_sw_20 {
            use clam::Metric;
            use algorithms_lib::Smith_Waterman::SmithWaterman;
            #[test]
            fn metric_sw_20to5_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test] 
            fn metric_sw_20to10_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_20to20_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_20to30_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_20to40_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_20to50_test(){
                let seq1_c: String = "GAAAGCCTATCGTCTGAGCG".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    
        mod test_sw_30 {
            use clam::Metric;
            use algorithms_lib::Smith_Waterman::SmithWaterman;
            #[test]
            fn metric_sw_30to5_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test] 
            fn metric_sw_30to10_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_30to20_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_30to30_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_30to40_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_30to50_test(){
                let seq1_c: String = "TGACCCACATTCATTCACTTATAGTATCTG".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    
        mod test_sw_40 {
            use clam::Metric;
            use algorithms_lib::Smith_Waterman::SmithWaterman;
            #[test]
            fn metric_sw_40to5_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "CAACC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test] 
            fn metric_sw_40to10_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TTGCTTTGAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_40to20_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_40to30_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_40to40_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
    
            #[test]
            fn metric_sw_40to50_test(){
                let seq1_c: String = "CCCTTAACCGAGTCTATGAACCTTATAAACAAGCCTTGCG".to_string();
                let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
                let seq1: String = seq1_c.clone();
                let seq2: String = seq2_c.clone();
                let metric: SmithWaterman = SmithWaterman{seq1, seq2};
                let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
            }
        }
    }
    
    mod test_sw_50 {
        use clam::Metric;
        use algorithms_lib::Smith_Waterman::SmithWaterman;
        #[test]
        fn metric_sw_50to5_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "CAACC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    
        #[test] 
        fn metric_sw_50to10_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TTGCTTTGAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    
        #[test]
        fn metric_sw_50to20_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "AAGGGACGCGTTGGAGTTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    
        #[test]
        fn metric_sw_50to30_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "GAGGAATCAAATCCTATGGCTATTGTCGGA".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman{seq1, seq2};
            let __score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    
        #[test]
        fn metric_sw_50to40_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "TCTGGTTCGATATTTCCCCGTTCTTCATCTACCGCCCTAC".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    
        #[test]
        fn metric_sw_50to50_test(){
            let seq1_c: String = "ATAATATCTAAAACAGCGATACTTGTATTGCTCGCCTGGGTACAATAGTT".to_string();
            let seq2_c: String = "ACTGATTTCCTTCTCACTCGATCAAATATTACCCGTGAGAGAATCGATAT".to_string();
            let seq1: String = seq1_c.clone();
            let seq2: String = seq2_c.clone();
            let metric: SmithWaterman = SmithWaterman{seq1, seq2};
            let _score:i8 = metric.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        }
    }
}