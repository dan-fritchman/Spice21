use std::sync::mpsc::{channel, Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;

pub trait Result {
    type Entry;
    fn new() -> Self;
    //  fn signals(&mut self, vars: &Variables<Entry>);
    fn push(&mut self, x: f64, vals: &Self::Entry);
    fn end(&mut self);
    //  fn len(&self) -> usize;
    //  fn get(&self, name: &str) -> SpResult<&Self::Entry>;
}

enum IoWriterMessage {
    STOP,
    DATA(Vec<f64>),
}
enum IoWriterResponse {
    OK,
    RESULT(Vec<Vec<f64>>),
}
/// Threaded File Writer
/// Temporarily somewhat abandoned
struct ThreadStreamWriter {
    tx: Sender<IoWriterMessage>,
    rx2: Receiver<IoWriterResponse>,
    thread: JoinHandle<()>,
}
impl ThreadStreamWriter {
    fn join(mut self) -> std::thread::Result<()> {
        self.thread.join()
    }
}
impl Result for ThreadStreamWriter {
    type Entry = Vec<f64>;
    fn new() -> Self {
        let (tx, rx) = channel::<IoWriterMessage>();
        let (tx2, rx2) = channel::<IoWriterResponse>();

        let thread = thread::spawn(move || {
            use serde::ser::{SerializeSeq, Serializer};
            use std::fs::File;

            let mut res: Vec<Vec<f64>> = vec![];

            let f = File::create("data.json").unwrap(); // FIXME: name
            let mut ser = serde_json::Serializer::new(f);
            let mut seq = ser.serialize_seq(None).unwrap();

            for msg in rx {
                match msg {
                    IoWriterMessage::DATA(d) => {
                        seq.serialize_element(&d).unwrap();
                        res.push(d);
                    }
                    IoWriterMessage::STOP => {
                        seq.end().unwrap();
                        tx2.send(IoWriterResponse::RESULT(res)).unwrap();
                        return;
                    }
                };
            }
        });
        Self { tx, rx2, thread }
    }
    fn push(&mut self, _x: f64, vals: &Vec<f64>) {
        self.tx.send(IoWriterMessage::DATA(vals.clone())).unwrap(); // FIXME: no copy
    }
    fn end(&mut self) {
        self.tx.send(IoWriterMessage::STOP).unwrap();
        for msg in self.rx2.iter() {
            match msg {
                IoWriterResponse::RESULT(_) => break,
                _ => continue,
            }
        }
    }
}
